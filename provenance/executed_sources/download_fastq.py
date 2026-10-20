#!/usr/bin/env python3
"""Retrieve and validate the original paired FASTQs; never use quality-reduced SRA Lite."""
import argparse
import concurrent.futures
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import urllib.parse

BASE = Path(__file__).resolve().parent
CURL = '/gpfs/ycga/work/coughlan/ps2327/ammo/envs/analysis/bin/curl'
DEFAULT_ROOT = Path('/vast/palmer/scratch/coughlan/ps2327/ammo_oct17_mapping/raw_fastq')


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(temp, path)


def discover(run):
    if not re.fullmatch(r'SRR\d+', run):
        raise ValueError(f'Invalid accession: {run}')
    url = 'https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?' + urllib.parse.urlencode(
        {'acc': run, 'accept-alternate-locations': 'yes'})
    result = subprocess.run([CURL, '--fail', '--silent', '--show-error', '--location',
                             '--retry', '5', '--retry-delay', '5', '--connect-timeout', '30',
                             '--max-time', '180', url], check=True, capture_output=True)
    body = json.loads(result.stdout)
    atomic_json(BASE / 'locator' / f'{run}.json', {'checked_utc': now(), 'url': url, 'body': body})
    bundles = [x for x in body.get('result', []) if x.get('bundle') == run]
    if len(bundles) != 1 or bundles[0].get('status') != 200:
        raise ValueError(f'{run}: locator did not return one successful bundle')
    files = [x for x in bundles[0].get('files', []) if x.get('type') == 'fastq']
    if len(files) != 2:
        raise ValueError(f'{run}: expected exactly two original FASTQ files, found {len(files)}')
    mates = {}
    for item in files:
        name = item['name']
        if Path(name).name != name or item.get('noqual'):
            raise ValueError(f'{run}: unsafe name or reduced-quality file')
        m = re.fullmatch(r'(.+?)(?:_R|_)([12])(?:_001)?\.(?:fq|fastq)\.gz', name)
        if not m or int(m[2]) in mates:
            raise ValueError(f'{run}: ambiguous paired FASTQ names: {name}')
        if int(item.get('size', 0)) <= 0 or not re.fullmatch('[0-9a-fA-F]{32}', item.get('md5', '')):
            raise ValueError(f'{run}: missing size or MD5 for {name}')
        links = [x['link'] for x in item.get('locations', []) if x.get('link', '').startswith('https://')]
        if not links:
            raise ValueError(f'{run}: no HTTPS original FASTQ location for {name}')
        mates[int(m[2])] = dict(item, mate=int(m[2]), stem=m[1], download_urls=links)
    if set(mates) != {1, 2} or mates[1]['stem'] != mates[2]['stem']:
        raise ValueError(f'{run}: FASTQ mate prefixes disagree')
    value = {'run': run, 'checked_utc': now(), 'source': url, 'files': [mates[1], mates[2]]}
    atomic_json(BASE / 'metadata' / f'{run}.json', value)
    return value


def validate(path, item):
    if path.stat().st_size != int(item['size']):
        raise ValueError(f'Size mismatch: {path}')
    digest = hashlib.md5()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            digest.update(block)
    if digest.hexdigest().lower() != item['md5'].lower():
        raise ValueError(f'MD5 mismatch: {path}')
    subprocess.run(['gzip', '-t', str(path)], check=True)
    return {'path': str(path), 'bytes': path.stat().st_size, 'md5': digest.hexdigest(),
            'gzip_validated': True}


def download(run, root):
    folder = root / run
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / '.download.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        meta_path = BASE / 'metadata' / f'{run}.json'
        metadata = json.loads(meta_path.read_text()) if meta_path.exists() else discover(run)
        verified = []
        for item in metadata['files']:
            final = folder / item['name']
            partial = folder / (item['name'] + '.partial')
            if final.exists():
                verified.append(dict(validate(final, item), mate=item['mate']))
                continue
            if partial.exists() and partial.stat().st_size > item['size']:
                raise ValueError(f'Oversized partial requires inspection: {partial}')
            # A full-size partial may be the complete file left by an interrupted validator.
            if not partial.exists() or partial.stat().st_size != item['size']:
                failure = None
                for url in item['download_urls']:
                    try:
                        subprocess.run([CURL, '--fail', '--location', '--show-error',
                                        '--continue-at', '-', '--retry', '8', '--retry-delay', '10',
                                        '--retry-all-errors', '--connect-timeout', '60',
                                        '--speed-time', '300', '--speed-limit', '1024',
                                        '--output', str(partial), url], check=True)
                        failure = None
                        break
                    except subprocess.CalledProcessError as error:
                        failure = error
                if failure:
                    raise failure
            info = validate(partial, item)
            os.replace(partial, final)
            info.update(path=str(final), mate=item['mate'])
            verified.append(info)
        marker = {'run': run, 'validated': True, 'verified_utc': now(),
                  'verification': 'NCBI original FASTQ size + MD5 + gzip CRC', 'files': verified}
        atomic_json(folder / 'VERIFIED.json', marker)
        atomic_json(BASE / 'verified' / f'{run}.json', marker)
        return marker


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--run')
    parser.add_argument('--index', type=int, help='Zero-based runs.txt index')
    parser.add_argument('--discover-all', action='store_true')
    parser.add_argument('--root', type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()
    runs = (BASE / 'runs.txt').read_text().splitlines()
    if args.discover_all:
        failures = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            future_runs = {pool.submit(discover, run): run for run in runs}
            for future in concurrent.futures.as_completed(future_runs):
                run = future_runs[future]
                try:
                    meta = future.result()
                    print(f'{run}: {sum(x["size"] for x in meta["files"])} bytes', flush=True)
                except Exception as error:
                    failures.append({'run': run, 'error': str(error)})
        atomic_json(BASE / 'DISCOVERY_STATUS.json', {'checked_utc': now(), 'runs': len(runs),
                                                   'failed': failures, 'complete': not failures})
        raise SystemExit(bool(failures))
    if (args.run is None) == (args.index is None):
        parser.error('Set exactly one of --run or --index')
    run = args.run if args.run else runs[args.index]
    if run not in runs:
        parser.error('Run is not in the saved project run list')
    try:
        print(json.dumps(download(run, args.root)), flush=True)
    except Exception as error:
        atomic_json(BASE / 'errors' / f'{run}.json', {'run': run, 'failed_utc': now(),
                                                    'error': repr(error)})
        raise


if __name__ == '__main__':
    main()
