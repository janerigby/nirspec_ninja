#!/usr/bin/env python3
"""Check, apply, or restore the JWST 3.0.0 FFT kernel-centering fix.

Run with the Python interpreter belonging to the JWST environment to patch.
Default action is read-only. Restart notebooks and pipeline workers after applying
or restoring. See https://github.com/spacetelescope/jwst/issues/11050.
"""
import argparse
import ast
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile

VERSION = '3.0.0'
ORIGINAL_CLASS_SHA = 'b265789baa2c1a03e3a1b02bd95cad29a4983c5909b1a4b4ae6ef3cbe37c07a1'
PATCHED_CLASS_SHA = '2b9b16d769d7dafe55cb1ac50aae9e32a053f8c75de5362586ca1d7cca80d8ec'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def class_text(source):
    tree = ast.parse(source)
    classes = [n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'NSClean']
    if len(classes) != 1:
        raise RuntimeError('Expected exactly one NSClean class. No changes made.')
    return ast.get_source_segment(source, classes[0]) + '\n'


def source_state(source):
    digest = sha(class_text(source).encode())
    if digest == ORIGINAL_CLASS_SHA:
        return 'original'
    if digest == PATCHED_CLASS_SHA:
        return 'patched'
    raise RuntimeError('NSClean source differs from both tested versions. No changes made. '
                       'Inspect the new source and the upstream fix before using this patch.')


def locate():
    version = importlib.metadata.version('jwst')
    if version != VERSION:
        raise RuntimeError(f'Installed JWST is {version}; this script supports only {VERSION}. '
                           'No changes made. Check whether the upgraded package includes the official fix.')
    spec = importlib.util.find_spec('jwst')
    if spec is None or spec.origin is None:
        raise RuntimeError('Cannot locate the installed JWST package.')
    path = Path(spec.origin).resolve().parent / 'clean_flicker_noise' / 'nsclean.py'
    if not path.is_file() or path.is_symlink():
        raise RuntimeError(f'Expected a regular source file: {path}')
    return path


def corrected(source):
    old = class_text(source)
    if source_state(source) != 'original':
        raise RuntimeError('Source is not the tested unmodified implementation.')
    if old.count('_mu = self.nx // 2 + 1') != 2 or old.count('_weight[self.ny // 2 + 1]') != 1:
        raise RuntimeError('Expected kernel-center expressions were not found.')
    new = old.replace('_mu = self.nx // 2 + 1', '_mu = self.nx // 2')
    new = new.replace('_weight[self.ny // 2 + 1]', '_weight[self.ny // 2]')
    result = source.replace(old.rstrip('\n'), new.rstrip('\n'), 1)
    if source_state(result) != 'patched':
        raise RuntimeError('Candidate patch does not match the tested correction.')
    return result


def verify_class(cls):
    """Independent numerical checks of BOTH kernels, for both orientations."""
    import numpy as np
    n = 128
    rng = np.random.default_rng(11050)
    mask = rng.random((n, n)) > .25
    mask[:, :4] = False
    mask[:, 20] = False
    impulse = np.zeros((n, n))
    impulse[40, 60] = 1
    lags = np.arange(-n // 2, n // 2)
    for detector in ('NRS1', 'NRS2'):
        obj = cls(detector, mask.copy(), weights_kernel_sigma=8)
        smooth = np.fft.irfft2(np.fft.rfft2(impulse) * obj.fgkern, s=impulse.shape)
        peak = tuple(map(int, np.unravel_index(smooth.argmax(), smooth.shape)))
        gaussian = np.exp(-lags**2 / (2 * obj.buffer_sigma**2)) / obj.buffer_sigma / np.sqrt(2 * np.pi)
        expected_smooth = sum(g * np.roll(impulse, int(lag), axis=0) for lag, g in zip(lags, gaussian))
        gaussian = np.exp(-lags**2 / (2 * 8**2)) / 8 / np.sqrt(2 * np.pi)
        expected_density = sum(g * np.roll(obj.mask.astype(float), int(lag), axis=1) for lag, g in zip(lags, gaussian))
        weights = obj.p_matrix[obj.mask]
        if peak != (40, 60) or not np.allclose(smooth, expected_smooth, rtol=0, atol=1e-6):
            raise RuntimeError(f'{detector}: smoothing kernel failed the zero-lag check.')
        if not np.all(np.isfinite(weights)) or not np.all(weights > 0):
            raise RuntimeError(f'{detector}: invalid background-sampling weights.')
        if not np.allclose(1 / weights, expected_density[obj.mask], rtol=0, atol=1e-6):
            raise RuntimeError(f'{detector}: weighting kernel failed direct-convolution check.')


def verify_candidate(source, path):
    # Execute the source being checked, not a previously imported/cached class.
    namespace = {'__name__': '_jwst_cfn_patch_check'}
    exec(compile(source, str(path), 'exec'), namespace)
    verify_class(namespace['NSClean'])


def verify_installed(path):
    # A new interpreter confirms that importing the installed package sees the fix.
    result = subprocess.run([sys.executable, '-B', str(Path(__file__).resolve()),
                             '--verify-installed', str(path)], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError('Fresh-process verification failed:\n' + result.stdout + result.stderr)
    print(result.stdout.strip())


def atomic_write(path, data, mode):
    fd, temporary = tempfile.mkstemp(prefix=path.name + '.tmp-', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, mode)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def backup_paths(path):
    return (path.with_name(path.name + '.cfn-lag-original.bak'),
            path.with_name(path.name + '.cfn-lag-backup.json'))


def check_backup(path, original, patched):
    backup, record = backup_paths(path)
    if backup.exists() or record.exists():
        if not backup.is_file() or not record.is_file():
            raise RuntimeError('Incomplete backup exists; inspect it before applying. No changes made.')
        metadata = json.loads(record.read_text())
        if backup.read_bytes() != original or metadata != {
            'jwst_version': VERSION, 'source_path': str(path),
            'original_sha256': sha(original), 'patched_sha256': sha(patched)}:
            raise RuntimeError('Existing backup does not match this installation. No changes made.')
        return
    # Never overwrite an existing backup.
    with backup.open('xb') as handle:
        handle.write(original)
    metadata = {'jwst_version': VERSION, 'source_path': str(path),
                'original_sha256': sha(original), 'patched_sha256': sha(patched)}
    with record.open('x') as handle:
        json.dump(metadata, handle, indent=2)
        handle.write('\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group()
    actions.add_argument('--check', action='store_true', help='Read-only status check (default).')
    actions.add_argument('--apply', action='store_true', help='Back up, apply the tested fix, and verify it.')
    actions.add_argument('--restore', action='store_true', help='Restore this script\'s original backup.')
    actions.add_argument('--verify-installed', metavar='PATH', help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    path = locate()
    original = path.read_bytes()
    state = source_state(original.decode('utf-8'))
    mode = stat.S_IMODE(path.stat().st_mode)
    if args.verify_installed:
        if path != Path(args.verify_installed).resolve() or state != 'patched':
            raise RuntimeError('Verification imported a different or unpatched installation.')
        from jwst.clean_flicker_noise.nsclean import NSClean
        verify_class(NSClean)
        print('PASS: installed weighting and smoothing kernels are centered (NRS1 and NRS2).')
        return
    print(f'Python: {sys.executable}\nJWST: {VERSION}\nSource: {path}\nStatus: {state}')
    if not args.apply and not args.restore:
        if state == 'patched':
            verify_candidate(original.decode('utf-8'), path)
            print('PASS: both corrected kernels verified. No files changed.')
        else:
            print('The tested one-pixel offsets are present. Use --apply to correct them.')
        return
    if args.restore:
        backup, record = backup_paths(path)
        if not backup.is_file() or not record.is_file():
            raise RuntimeError('No complete backup from this script is available. No changes made.')
        saved = backup.read_bytes()
        metadata = json.loads(record.read_text())
        if (metadata.get('jwst_version') != VERSION or metadata.get('source_path') != str(path)
                or sha(saved) != metadata.get('original_sha256')
                or source_state(saved.decode('utf-8')) != 'original'):
            raise RuntimeError('Backup verification failed. No changes made.')
        if state == 'original' and original == saved:
            print('Already restored. No files changed.')
            return
        if state != 'patched' or sha(original) != metadata.get('patched_sha256'):
            raise RuntimeError('Installed file has other changes; refusing to overwrite them.')
        if path.read_bytes() != original:
            raise RuntimeError('Source changed during the check. No changes made.')
        atomic_write(path, saved, mode)
        print('Restored the original JWST 3.0.0 source. Backup retained.')
    elif state == 'patched':
        verify_installed(path)
        print('The tested fix is already present. No files changed.')
        return
    else:
        patched = corrected(original.decode('utf-8')).encode('utf-8')
        verify_candidate(patched.decode('utf-8'), path)
        check_backup(path, original, patched)
        if path.read_bytes() != original:
            raise RuntimeError('Source changed during the check. No source changes made.')
        atomic_write(path, patched, mode)
        try:
            verify_installed(path)
        except Exception:
            # Roll back only if nobody has changed the file since our write.
            if path.read_bytes() == patched:
                atomic_write(path, original, mode)
                print('Verification failed; the original source was restored.', file=sys.stderr)
            raise
        print(f'Applied the three-line fix. Backup: {backup_paths(path)[0]}')
    print('Restart running notebooks, Python sessions, and pipeline workers before using this environment.')


if __name__ == '__main__':
    try:
        # Keep incidental matplotlib cache creation outside the installation.
        with tempfile.TemporaryDirectory(prefix='jwst-cfn-patch-') as cache:
            os.environ.setdefault('MPLCONFIGDIR', cache)
            main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
