"""Process scheduling and CRDS preparation, separate from stage execution.

Importing this module does not configure CRDS or start workers.
"""
import os
from pathlib import Path
import subprocess
import sys
from multiprocessing.pool import Pool

_crds_preflight_context = None

def preflight_crds(input_files, stage_name):
    """Cache CRDS mappings and input-specific best references before a pool."""
    global _crds_preflight_context
    if _crds_preflight_context is None:
        import crds
        from crds.client import api as crds_api

        context = crds.get_default_context('jwst')
        print(f'CRDS preflight: caching mapping closure for {context}')
        crds_api.dump_mappings(context)

        # Spawned workers inherit the automatically selected literal context
        # and do not repeat context resolution or race to cache mappings.
        os.environ['CRDS_CONTEXT'] = context
        _crds_preflight_context = context
        print(f'CRDS mapping preflight complete: {context}')

    context = _crds_preflight_context
    fits_files = sorted({
        os.path.abspath(filename) for filename in input_files
        if str(filename).lower().endswith(('.fits', '.fit'))
    })
    if not fits_files:
        print(f'CRDS bestrefs preflight: no FITS inputs for {stage_name}')
        return context

    crds_executable = os.path.join(os.path.dirname(sys.executable), 'crds')
    command = [
        crds_executable, 'bestrefs',
        '--new-context', context,
        '--sync-mappings', '1',
        '--sync-references', '1',
        '--files', *fits_files,
    ]
    print(
        f'CRDS bestrefs preflight: checking {len(fits_files)} '
        f'{stage_name} input file(s).\n'
        'CRDS will now download any missing reference files. Some NIRSpec '
        'references are several GB, so this may take a while and may appear '
        'inactive. The worker pool will start only after all downloads finish.',
        flush=True,
    )
    subprocess.run(command, check=True)
    print(
        f'CRDS bestrefs preflight complete for {stage_name}; starting workers.',
        flush=True,
    )
    return context



def run_pipeline_tasks(worker, tasks, cores2use):
    """Run prepared ``(input_file, parameters)`` tasks sequentially or in a pool."""
    if not tasks:
        return []
    if not isinstance(cores2use, int) or isinstance(cores2use, bool) or cores2use < 1:
        raise ValueError(f'cores2use must be a positive integer; got {cores2use!r}')
    if cores2use == 1:
        return [worker(*task) for task in tasks]

    worker_count = min(cores2use, len(tasks))
    pool = Pool(worker_count)
    try:
        results = pool.starmap(worker, tasks)
    except BaseException:
        pool.terminate()
        pool.join()
        raise
    else:
        pool.close()
        pool.join()
        return results



def predownload_crds(myfiles, context="jwst_1584.pmap"):
    # Download all needed JWST CRDS files *before* calling the pipeline.  Prevents collisions when multiprocessing
    subprocess.run(["crds", "bestrefs","--new-context", context,
        "--files", *map(str, myfiles), "--sync-references=1",], check=True)



def max_filesize(indir):
    folder_path = Path(indir)
    # Find the largest filesize for a bunch of files in a directory (ignore directories)
    files = [f for f in folder_path.rglob("*") if f.is_file()]
    if files:
        largest_file = max(files, key=lambda f: f.stat().st_size)
        return(largest_file.stat().st_size)
    else: raise Exception("No files found")

