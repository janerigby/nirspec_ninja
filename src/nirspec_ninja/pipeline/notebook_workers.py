"""Existing fixed-slit notebook stage workers, migrated without changing settings."""
from nirspec_ninja.pipeline.execution import close_pipeline_result as _close_pipeline_result

from os.path import basename
import traceback

def simple_task_to_multiprocess(arg1, arg2, arg3):
    return arg1 + 1


def _record_failure(directory, filename, stage): 
    from pathlib import Path 
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True) 
    path = directory / (filename + '_pipecrash.txt') 
    path.write_text(traceback.format_exc()) 
    print(f'\n *** {stage} failed ***') 
    traceback.print_exc() 
    print('Full traceback saved to:', path) 

    
def run_jwst_det1(uncal_file, output_dir, paramdict, save_rateints=False):
    """Run Detector1 with cleaning skipped; return RATE path or None.

    Explicit step settings take precedence over fallback defaults. This notebook
    workflow requires ramp fitting and saved RATE output, with cleaning in Spec2.
    """
    from copy import deepcopy 
    from pathlib import Path 
    from time import perf_counter
    started = perf_counter() 
    source = Path(uncal_file).expanduser().resolve()
    outdir = Path(output_dir).expanduser().resolve()
    success = False 
    try: 
        from jwst.pipeline.calwebb_detector1 import Detector1Pipeline 
        outdir.mkdir(parents=True, exist_ok=True) 
        params = deepcopy(paramdict) 
        if 'clean_flicker_noise' in Detector1Pipeline.step_defs: 
            cleaning = params.setdefault('clean_flicker_noise', {})
            cleaning.setdefault('skip', True)
            if not cleaning['skip']:
                raise ValueError('This workflow requires Detector1 clean_flicker_noise.skip=True; cleaning runs in Spec2')
        # Ramp fitting is required to produce the RATE files used by later stages.
        params.setdefault('ramp_fit', {}).setdefault('skip', False)
        if params.get('ramp_fit', {}).get('skip', False): 
            raise ValueError('This workflow requires ramp_fit to produce a RATE file') 
        from nirspec_ninja.pipeline.detector1_outputs import configure_rateints_saving
        # Save the uncleaned RATE files required by the mask and Spec2 stages.
        config, config_file = Detector1Pipeline.build_config(
            str(source), output_dir=str(outdir), steps=params, save_results=True,
        )
        pipeline = Detector1Pipeline.from_config_section(config, config_file=config_file)
        configure_rateints_saving(pipeline, save_rateints=save_rateints)
        result = pipeline.run(str(source))
        _close_pipeline_result(result) 
        rate_path = outdir / (source.stem.removesuffix('_uncal') + '_rate.fits') 
        if not rate_path.is_file(): 
            raise FileNotFoundError(f'Detector1 did not produce {rate_path}')
        success = True 
        #print('Detector1 RATE ready:', rate_path) 
        return str(rate_path) 
    except Exception: 
        _record_failure(outdir, source.stem, 'Detector1') 
        return None 
    finally:
        print(f'Exposure {source.name} | {perf_counter() - started:.1f} s | ' 
              f'{"SUCCESS" if success else "FAILURE"}', flush=True) 


def make_rate_custom_mask(rate_file, input_dir, paramdict=None):
    (
        'Mask an existing uncleaned RATE in its original mode; return mask path or None.\n' 
        '\n' 
        '    input_dir locates the MSA metadata beside the original UNCAL files.\n'
    )
    from copy import deepcopy
    from pathlib import Path 
    source = Path(rate_file).expanduser().resolve() 
    mask_dir = source.parent / 'custom_masks' 
    try:
        from jwst import datamodels
        if not source.name.endswith('_rate.fits'): 
            raise ValueError(f'Expected a *_rate.fits input: {source}') 
        with datamodels.open(source) as rate: 
            if str(rate.meta.cal_step.clean_flicker_noise).upper() == 'COMPLETE':
                raise ValueError(f'RATE already flicker-cleaned: {source}. Regenerate with Detector1 cleaning skipped.') 
        from nirspec_ninja.masks.custom_nirspec_mask_2 import make_custom_NIRSpec_mask 
        mask_dir.mkdir(parents=True, exist_ok=True) 
        mask_path = mask_dir / (source.stem.removesuffix('_rate') + '_mask.fits') 
        params = dict(mask_MSA_by_wcs=True, mask_FS_by_wcs=True,
                      mask_by_illumination=True, mask_FS_strip=False, mask_FSprism=False,
                      mask_illumination_nsigma=5.0, mask_grating_extrapolate=True)
        if paramdict is not None:
            params.update(deepcopy(paramdict))
        make_custom_NIRSpec_mask(str(source), str(mask_path), **params,
                                overwrite=True, input_dir=str(Path(input_dir).expanduser().resolve()))
        if not mask_path.is_file(): 
            raise FileNotFoundError(f'Mask generation did not produce {mask_path}') 
        print(mask_path) 
        return str(mask_path) 
    except Exception:
        _record_failure(mask_dir, source.stem + '_mask', 'RATE mask generation') 
        return None 


def run_jwst_det1_custom_mask(uncal_file, output_dir, paramdict):
    """Optional combined worker; notebook uses the two separate stages instead.""" 
    from pathlib import Path 
    rate_file = run_jwst_det1(uncal_file, output_dir, paramdict) 
    if rate_file is None:   return None 
    mask_file = make_rate_custom_mask(rate_file, Path(uncal_file).expanduser().resolve().parent) 
    return rate_file if mask_file is not None else None 


def run_extract(L3_destripedir, L3_extractdir, clip_nan=5):
    """Return a compact directory report, including per-file trace failures."""
    from pathlib import Path
    from nirspec_ninja.fixed_slit.extraction import wrap_median_combine_level3_nirspecFS
    report = {'directory': str(Path(L3_destripedir).resolve()),
              'output_directory': str(Path(L3_extractdir).resolve()),
              'status': 'ok', 'files': [], 'reason': None}
    try:
        report['files'] = wrap_median_combine_level3_nirspecFS(
            L3_destripedir, L3_extractdir, clip_nan=clip_nan, return_report=True)
        if any(item['status'] == 'failed' for item in report['files']):
            report['status'] = 'failed'
        elif not report['files']:
            report['status'] = 'empty'
            report['reason'] = 'no eligible S2D/CAL files (S200B1 is excluded)'
            print(f"EXTRACTION EMPTY | {report['directory']} | {report['reason']}", flush=True)
    except Exception as error:
        report['status'] = 'failed'
        report['reason'] = f'{type(error).__name__}: {error}'
        print(f"EXTRACTION FAILED | {report['directory']} | {report['reason']}", flush=True)
    return report

def summarize_extract(results):
    """Print one short summary; workers already identify affected paths."""
    files = [item for report in results for item in report['files']]
    written = sum(item['status'] == 'ok' for item in files)
    problems = []
    counts = (
        ('extraction errors', sum(item['status'] == 'failed' for item in files)
         + sum(report['status'] == 'failed' and report['reason'] is not None
               for report in results)),
        ('trace failures', sum(item.get('trace', {}).get('status') == 'failed'
                               for item in files)),
        ('empty directories', sum(report['status'] == 'empty' for report in results)),
    )
    for label, count in counts:
        if count:
            problems.append(f'{count} {label}')
    suffix = '; ' + ', '.join(problems) if problems else '; no failures'
    print(f'Extraction: {written} CSVs written{suffix}.')


def run_jwst_JR_custom_extract(L3_destripedir, L3_extractdir, clip_nan=5):
    from nirspec_ninja.fixed_slit.extraction import wrap_median_combine_level3_nirspecFS
    spec1 = wrap_median_combine_level3_nirspecFS(L3_destripedir, L3_extractdir, clip_nan)
    return(spec1)

        
def run_jwst_spec3(asn_file, output_dir, paramdict=None):
    import warnings
    warnings.filterwarnings("ignore", message=r"Input association file contains path information;.*", \
        category=UserWarning, module=r"jwst\.associations\.association",)
    from copy import deepcopy
    from jwst.pipeline import Spec3Pipeline
    from nirspec_ninja.pipeline.execution import apply_pipeline_parameters
    spec3 = Spec3Pipeline()
    spec3.outlier_detection.skip = False
    if paramdict is not None:
        apply_pipeline_parameters(spec3, deepcopy(paramdict))
    spec3.output_dir = output_dir
    spec3.save_results = True
    print("Starting Spec3:", output_dir, flush=True)
    _close_pipeline_result(spec3.run(asn_file))


def save_cleaned_rate_in_l2a(clean_step, source, exposure_type, opmode):
    """Save the Stage 2 cleaning result before subsequent calibration steps."""
    from pathlib import Path
    source = Path(source)
    path = source.with_name(source.stem.removesuffix('_rate') + '_rate_after_cleaning.fits')
    original_run = clean_step.run

    def run_and_save(*args, **kwargs):
        result = original_run(*args, **kwargs)
        status = str(result.meta.cal_step.clean_flicker_noise).upper()
        if status in ('NONE', 'SKIPPED'):
            status = str(getattr(result.meta.cal_step, 'nsclean', None)).upper()
        if status != 'COMPLETE':
            raise RuntimeError(f'Cleaning did not complete for {source}: {status}')
        # The FS view is used downstream; the diagnostic retains source metadata.
        with result.copy() as diagnostic:
            diagnostic.meta.exposure.type = exposure_type
            if opmode is None:
                # Do not carry the temporary FS mode into an originally absent field.
                if hasattr(diagnostic.meta.instrument, 'opmode'):
                    del diagnostic.meta.instrument.opmode
            else:
                diagnostic.meta.instrument.opmode = opmode
            diagnostic.save(str(path))
        return result

    clean_step.run = run_and_save


def run_jwst_spec2(infile, output_dir, paramdict, save_results=True):
    """Clean with the matching mask, then calibrate an in-memory fixed-slit view.

    Explicit cleaning settings take precedence over fallback defaults. Cleaning
    must be enabled. The source RATE stays unchanged; its cleaned diagnostic is
    saved independently of pipeline-level save_results. Calibrated products are
    saved by default for Spec3; save_results=False is available for diagnostic
    callers. The combined NSClean patch is installed and verified in this
    worker with its default weight cap enabled before cleaning. Return success bool.
    """
    from copy import deepcopy 
    from pathlib import Path
    source = Path(infile).expanduser().resolve() 
    outdir = Path(output_dir).expanduser().resolve() 
    try:
        from jwst import datamodels 
        from jwst.pipeline import Spec2Pipeline 
        mask_path = source.parent / 'custom_masks' / (source.stem.removesuffix('_rate') + '_mask.fits') 
        if not mask_path.is_file():
            raise FileNotFoundError(f'Spec2 requires {mask_path}. Run the notebook mask stage first.') 
        params = deepcopy(paramdict) 
        cleaning_name = 'clean_flicker_noise' if 'clean_flicker_noise' in Spec2Pipeline.step_defs else 'nsclean' 
        other_name = 'nsclean' if cleaning_name == 'clean_flicker_noise' else 'clean_flicker_noise' 
        fallback = params.pop(other_name, {}) 
        cleaning = params.setdefault(cleaning_name, fallback) 
        cleaning.setdefault('skip', False)
        if cleaning['skip']:
            raise ValueError('This workflow requires Spec2 clean_flicker_noise.skip=False')
        cleaning['user_mask'] = str(mask_path)
        cleaning.setdefault('save_results', False)
        cleaning.setdefault('save_noise', False)
        cleaning.setdefault('save_background', False)
        cleaning.setdefault('fit_method', 'median') 
        cleaning.setdefault('background_method', None) 
        cleaning.setdefault('fit_by_channel', False) 
        cleaning.setdefault('autoparam', False) 
        with datamodels.open(source) as rate: 
            if str(rate.meta.cal_step.clean_flicker_noise).upper() == 'COMPLETE': 
                raise ValueError(f'RATE already flicker-cleaned: {source}. Regenerate with Detector1 cleaning skipped.') 
            exposure_type = rate.meta.exposure.type
            opmode = getattr(rate.meta.instrument, 'opmode', None)
            rate.meta.exposure.type = 'NRS_FIXEDSLIT'
            rate.meta.instrument.opmode = 'FIXEDSLIT'
            outdir.mkdir(parents=True, exist_ok=True)
            config, config_file = Spec2Pipeline.build_config(
                rate, input_dir=str(source.parent), output_dir=str(outdir),
                steps=params, save_results=save_results,
            )
            pipeline = Spec2Pipeline.from_config_section(config, config_file=config_file)
            from nirspec_ninja.pipeline.patches.nirspec_flicker_patch import (
                install_patch, save_patch_record, DEFAULT_MAX_WEIGHT)
            patch_status = install_patch(max_weight=DEFAULT_MAX_WEIGHT)
            getattr(pipeline, cleaning_name).flicker_patch_status = patch_status
            save_patch_record(mask_path, patch_status)
            save_cleaned_rate_in_l2a(getattr(pipeline, cleaning_name), source,
                                     exposure_type, opmode)
            result = pipeline.run(rate)
            _close_pipeline_result(result)
            from nirspec_ninja.diagnostics.flicker_noise_diagnostics import write_noise_report
            cleaned_rate = source.with_name(
                source.stem.removesuffix('_rate') + '_rate_after_cleaning.fits')
            write_noise_report(source, cleaned_rate, mask_path,
                report_file=source.parent / 'flicker_noise_rms.txt',
                fit_method=cleaning['fit_method'])
        return True 
    except Exception:
        _record_failure(outdir, source.stem + '_spec2', 'Spec2') 
        return False
