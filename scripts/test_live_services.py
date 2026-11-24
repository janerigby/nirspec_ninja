#!/usr/bin/env python3
"""Test live WISE retrieval and one JWST product downloaded from AWS and MAST.

Run in your JWST reduction Python environment. Supply a NIRSpec exposure with
pointing metadata. No pipeline reduction or notebook modification is performed.
"""
import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path
import subprocess
import sys
import tempfile


def checksum(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def command(args, log, timeout):
    with log.open('w') as handle:
        result = subprocess.run(args, stdout=handle, stderr=subprocess.STDOUT,
                                timeout=timeout)
    if result.returncode:
        raise RuntimeError(f'Command failed (exit {result.returncode}); see {log}')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('exposure', type=Path,
                        help='Local NIRSpec exposure FITS file with pointing metadata')
    parser.add_argument('--archive-file',
                        help='Public JWST FITS basename; default: RATE corresponding to exposure')
    parser.add_argument('--output', type=Path,
                        help='New output directory; default: new directory under AI_experiments/Temp')
    parser.add_argument('--max-download-mb', type=float, default=512,
                        help='Maximum size per archive copy in decimal MB (default 512)')
    parser.add_argument('--timeout', type=int, default=900,
                        help='Time limit for each download/chart command, seconds (default 900)')
    args = parser.parse_args(argv)
    exposure = args.exposure.expanduser().resolve()
    if not exposure.is_file():
        parser.error(f'Exposure does not exist: {exposure}')
    if args.timeout <= 0 or args.max_download_mb <= 0:
        parser.error('Timeout and download limit must be positive')
    filename = args.archive_file or exposure.name.replace('_uncal.fits', '_rate.fits')
    if Path(filename).name != filename or not filename.endswith('.fits'):
        parser.error('--archive-file must be a FITS basename, not a path')
    if args.output:
        output = args.output.expanduser().resolve()
        output.mkdir(parents=True, exist_ok=False)
    else:
        parent = Path.home() / 'AI_experiments/Temp'
        parent.mkdir(parents=True, exist_ok=True)
        output = Path(tempfile.mkdtemp(prefix='nirspec_live_services_', dir=parent))
    print(f'Outputs: {output}', flush=True)
    summary = {'started': datetime.now().astimezone().isoformat(),
               'exposure': str(exposure), 'archive_file': filename, 'tests': {}}

    def save():
        (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')

    save()
    print('Testing WISE retrieval and finder chart...', flush=True)
    try:
        from astropy.io import fits
        from nirspec_ninja.finder.nirspec_finder import geometry
        geometry(exposure)  # Check pointing metadata before contacting IRSA.
        chart = output / 'wise_finder.png'
        command([sys.executable, '-m', 'nirspec_ninja.finder.nirspec_finder',
                 str(exposure), '--output', str(chart), '--cache',
                 str(output / 'wise_cache'), '--regions'],
                output / 'wise.log', args.timeout)
        record = json.loads(chart.with_suffix('.json').read_text())
        if not chart.is_file() or not chart.with_suffix('.reg').is_file():
            raise RuntimeError('Finder command did not produce PNG and region files')
        summary['tests']['WISE'] = {'status': 'pass', 'chart': str(chart),
                                   'image_source': record['image']}
        print('WISE: PASS', flush=True)
    except Exception as error:
        summary['tests']['WISE'] = {'status': 'fail', 'error': str(error)}
        print(f'WISE: FAIL — {error}', flush=True)
    save()

    inventory = output / 'filenames.txt'
    inventory.write_text(filename + '\n')
    hashes = {}
    for source in ('AWS', 'MAST'):
        print(f'Testing {source} download...', flush=True)
        try:
            import requests
            from astropy.io import fits
            from nirspec_ninja.archive import MAST_download_rawfilelist as downloader
            downloader.SOURCE_POLICY = source
            url = downloader.urls_for(filename)[0][1]
            with requests.Session() as session:
                remote = downloader.inspect_remote(session, source, url, filename)
            size_mb = remote.size / 1e6
            if size_mb > args.max_download_mb:
                raise RuntimeError(f'File is {size_mb:.1f} MB, above the '
                                   f'{args.max_download_mb:g} MB limit. Choose a smaller '
                                   '--archive-file or explicitly raise --max-download-mb.')
            print(f'{source}: {size_mb:.1f} MB; downloading one file', flush=True)
            destination = output / source.lower()
            command([sys.executable, '-m', 'nirspec_ninja.archive.MAST_download_rawfilelist',
                     str(inventory), '--source', source, '--workers', '1',
                     '--output-dir', str(destination)], output / f'{source.lower()}.log',
                    args.timeout)
            path = destination / filename
            if path.stat().st_size != remote.size:
                raise RuntimeError('Downloaded size differs from remote size')
            with fits.open(path) as hdus:
                hdus.verify('exception')
                if hdus[0].header.get('INSTRUME', '').upper() != 'NIRSPEC':
                    raise RuntimeError('Downloaded file does not identify NIRSPEC')
            digest = checksum(path)
            hashes[source] = digest
            summary['tests'][source] = {'status': 'pass', 'file': str(path),
                                        'bytes': path.stat().st_size, 'sha256': digest}
            print(f'{source}: PASS', flush=True)
        except Exception as error:
            summary['tests'][source] = {'status': 'fail', 'error': str(error)}
            print(f'{source}: FAIL — {error}', flush=True)
        save()
    if len(hashes) == 2:
        equal = hashes['AWS'] == hashes['MAST']
        summary['tests']['AWS_vs_MAST'] = {'status': 'pass' if equal else 'fail',
                                          'byte_identical': equal}
        print(f'AWS/MAST byte comparison: {"PASS" if equal else "FAIL"}', flush=True)
    summary['finished'] = datetime.now().astimezone().isoformat()
    save()
    print(f'Summary: {output / "summary.json"}', flush=True)
    return int(any(item['status'] != 'pass' for item in summary['tests'].values()))


if __name__ == '__main__':
    raise SystemExit(main())
