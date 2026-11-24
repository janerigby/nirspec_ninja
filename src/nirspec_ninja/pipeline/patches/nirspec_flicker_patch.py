"""Verified, process-local JWST 3.0.0 NSClean patch; no site-package writes.

Centers both kernels, bounds inverse-density weights, and rejects saturated
ramp groups. Call install_patch in every worker before running cleaning.
"""
import hashlib
import inspect
import json
from pathlib import Path

import jwst
import numpy as np
from stdatamodels.jwst.datamodels import dqflags
from jwst.clean_flicker_noise import nsclean, clean_flicker_noise as cfn

PATCH_REVISION = 'nirspec-flicker-dq-cap-v1'
DEFAULT_MAX_WEIGHT = 4.0
PIXEL_EXCLUSIONS = ('DO_NOT_USE', 'NONLINEAR', 'NO_LIN_CORR', 'NO_SAT_CHECK')
_ORIGINAL_HASH = 'b265789baa2c1a03e3a1b02bd95cad29a4983c5909b1a4b4ae6ef3cbe37c07a1'
_CENTERED_HASH = '2b9b16d769d7dafe55cb1ac50aae9e32a053f8c75de5362586ca1d7cca80d8ec'
_BASE_CLASS = nsclean.NSClean
_BASE_GROUP_MASK = cfn._mask_unusable
_CENTERED_CLASS = None


def validate_cap(max_weight):
    if isinstance(max_weight, (bool, np.bool_)):
        raise ValueError('NSClean max_weight must be finite and >= 1')
    value = float(max_weight)
    if not np.isfinite(value) or value < 1:
        raise ValueError('NSClean max_weight must be finite and >= 1')
    return value


def _centered_class():
    global _CENTERED_CLASS
    if _CENTERED_CLASS is not None:
        return _CENTERED_CLASS
    if jwst.__version__ != '3.0.0':
        raise RuntimeError(f'{PATCH_REVISION} supports tested JWST 3.0.0 only; '
                           f'found {jwst.__version__}. Review upstream before upgrading.')
    source = inspect.getsource(_BASE_CLASS).rstrip() + '\n'
    digest = hashlib.sha256(source.encode()).hexdigest()
    if digest not in (_ORIGINAL_HASH, _CENTERED_HASH):
        raise RuntimeError('NSClean source differs from both tested JWST 3.0.0 '
                           'versions; refusing to apply an unverified patch.')
    if digest == _ORIGINAL_HASH:
        source = source.replace('_mu = self.nx // 2 + 1', '_mu = self.nx // 2')
        source = source.replace('_weight[self.ny // 2 + 1]', '_weight[self.ny // 2]')
    if hashlib.sha256(source.encode()).hexdigest() != _CENTERED_HASH:
        raise RuntimeError('Centered NSClean source verification failed')
    namespace = vars(nsclean).copy()
    exec(compile(source, '<verified-centered-NSClean>', 'exec'), namespace)
    _CENTERED_CLASS = namespace['NSClean']
    return _CENTERED_CLASS


def _mask_unusable_with_saturation(mask, dq):
    _BASE_GROUP_MASK(mask, dq)
    mask[(dq & dqflags.group['SATURATED']) != 0] = False


def verify_patch(max_weight):
    """Check active bindings, both kernel orientations, cap, and group DQ."""
    cap = validate_cap(max_weight)
    cls = cfn.NSClean
    if (cls is not nsclean.NSClean or getattr(cls, 'patch_revision', None) != PATCH_REVISION
            or cls.configured_max_weight != cap
            or cfn._mask_unusable is not _mask_unusable_with_saturation):
        raise RuntimeError('NSClean runtime patch is not active in this process')
    n = 128
    mask = np.ones((n, n), bool)
    mask[20:100, 30:90] = False
    mask[59, 60] = True  # isolated admitted sample
    impulse = np.zeros((n, n)); impulse[40, 60] = 1
    lags = np.arange(-n // 2, n // 2)
    for detector in ('NRS1', 'NRS2'):
        obj = cls(detector, mask.copy(), weights_kernel_sigma=8)
        gaussian = np.exp(-lags**2 / (2 * 8**2)) / 8 / np.sqrt(2 * np.pi)
        density = sum(g * np.roll(obj.mask.astype(float), int(lag), axis=1)
                      for lag, g in zip(lags, gaussian))
        expected = 1 / np.maximum(density[obj.mask], 1 / cap)
        actual = obj.p_matrix[obj.mask]
        if not np.allclose(actual, expected, rtol=2e-6, atol=2e-6):
            raise RuntimeError(f'{detector}: centered density/cap check failed')
        if not np.all((actual > 0) & (actual <= cap)):
            raise RuntimeError(f'{detector}: invalid or unbounded weights')
        smooth = np.fft.irfft2(np.fft.rfft2(impulse) * obj.fgkern, s=impulse.shape)
        gaussian = np.exp(-lags**2 / (2 * obj.buffer_sigma**2)) / obj.buffer_sigma / np.sqrt(2 * np.pi)
        expected_smooth = sum(g * np.roll(impulse, int(lag), axis=0)
                              for lag, g in zip(lags, gaussian))
        if not np.allclose(smooth, expected_smooth, atol=1e-6, rtol=0):
            raise RuntimeError(f'{detector}: smoothing kernel check failed')
    flags = np.array([0, dqflags.group['DO_NOT_USE'], dqflags.group['SATURATED'],
                      dqflags.group['JUMP_DET']], dtype=np.uint8)
    usable = np.ones(4, bool); cfn._mask_unusable(usable, flags)
    if not np.array_equal(usable, [True, False, False, False]):
        raise RuntimeError('Group-DQ exclusion check failed')
    return dict(patch_revision=PATCH_REVISION, jwst_version=jwst.__version__,
                max_weight=cap, kernel_centers_verified=True,
                saturated_groups_excluded=True, runtime_verified=True,
                pixel_exclusions=list(PIXEL_EXCLUSIONS))


def install_patch(max_weight=DEFAULT_MAX_WEIGHT):
    """Apply/reconfigure and numerically verify the patch in the calling process."""
    cap = validate_cap(max_weight)
    base = _centered_class()
    class BoundedNSClean(base):
        patch_revision = PATCH_REVISION
        configured_max_weight = cap
        def __init__(self, *args, max_weight=cap, **kwargs):
            self.max_weight = validate_cap(max_weight)
            super().__init__(*args, **kwargs)
            np.minimum(self.p_matrix, self.max_weight, out=self.p_matrix)
    previous = nsclean.NSClean, cfn.NSClean, cfn._mask_unusable
    cfn.NSClean = nsclean.NSClean = BoundedNSClean
    cfn._mask_unusable = _mask_unusable_with_saturation
    try:
        return verify_patch(cap)
    except Exception:
        nsclean.NSClean, cfn.NSClean, cfn._mask_unusable = previous
        raise


def save_patch_record(mask_file, status):
    path = Path(mask_file).with_suffix('.flicker_patch.json')
    path.write_text(json.dumps(status, indent=2) + '\n')
    return path


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Apply and verify the patch in this process only')
    parser.add_argument('--max-weight', type=float, default=DEFAULT_MAX_WEIGHT)
    args = parser.parse_args()
    print(json.dumps(install_patch(args.max_weight), indent=2))
