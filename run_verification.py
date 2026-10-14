"""Reproduce the public checks without the manuscript or private files.

Default: verify payload hashes/notebooks and run current core Python checks.
--mathematica: additionally run all 5 standalone WL programs in separate
temporary directories and fresh kernels. Logs and receipt go in --output.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent


def normalized_bytes(path):
    return path.read_bytes().replace(b'\r\n', b'\n')


def verify_payload():
    entries = [line for line in (ROOT/'MANIFEST.sha256').read_text().splitlines()
               if line and not line.startswith('#')]
    for line in entries:
        expected, name = line.split('  ', 1)
        actual = hashlib.sha256(normalized_bytes(ROOT/name)).hexdigest()
        if actual.lower() != expected.lower():
            raise ValueError('Payload hash mismatch: '+name)
    for wl in sorted((ROOT/'mathematica').glob('*.wl')):
        nb = wl.with_suffix('.nb').read_text(encoding='utf-8')
        matches = re.findall(r'Cell\[("(?:\\.|[^"\\])*"),"Program"\]', nb)
        if len(matches) != 1:
            raise ValueError('Expected one readable source cell: '+wl.name)
        payload = json.loads(matches[0])
        if payload.replace('\r\n','\n').rstrip() != wl.read_text(encoding='utf-8').rstrip():
            raise ValueError('Notebook/WL mismatch: '+wl.name)
    print(f'Payload hashes: {len(entries)} matched; 5 notebook payloads matched.')


def run_one(command, cwd, output, timeout):
    started = time.monotonic()
    completed = subprocess.run(command, cwd=cwd, capture_output=True,
                               text=True, encoding='utf-8', errors='replace', timeout=timeout)
    log = completed.stdout+completed.stderr
    output.write_text(log, encoding='utf-8')
    return completed.returncode, log, round(time.monotonic()-started,3)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mathematica', action='store_true')
    parser.add_argument('--wolframscript', default='wolframscript')
    parser.add_argument('--output', type=Path, default=ROOT/'verification-output')
    parser.add_argument('--timeout', type=int, default=1800)
    args = parser.parse_args()
    verify_payload()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True,exist_ok=True)
    results = []
    # Also demonstrate that the Python calculation has no repository inputs.
    with tempfile.TemporaryDirectory(prefix='nrh-core-', dir=args.output) as directory:
        copied = Path(directory)/'verify_core.py'
        shutil.copyfile(ROOT/'checks'/copied.name,copied)
        code,log,seconds = run_one([sys.executable,str(copied)], directory,
                                   args.output/'core-python.txt', args.timeout)
    passed = log.count('[PASS]'); failed = log.count('[FAIL]')
    results.append(dict(program='verify_core.py',exit=code,passed=passed,
                        failed=failed,seconds=seconds,status='passed' if code==0 and passed>0 and failed==0 else 'failed'))
    if args.mathematica:
        executable = shutil.which(args.wolframscript)
        if executable is None:
            raise FileNotFoundError('wolframscript not found; provide --wolframscript PATH')
        for wl in sorted((ROOT/'mathematica').glob('*.wl')):
            with tempfile.TemporaryDirectory(prefix='nrh-wl-', dir=args.output) as directory:
                copied = Path(directory)/wl.name
                shutil.copyfile(wl,copied)
                code,log,seconds = run_one([executable,'-file',str(copied)], directory,
                                           args.output/(wl.stem+'.txt'),args.timeout)
            summaries = re.findall(r'(\d+)\s*/\s*(\d+)\s+checks passed',log)
            passed,total = map(int,summaries[-1]) if summaries else (0,0)
            diagnostic = bool(re.search(r'\b\w+::\w+\s*:',log))
            ok = code==0 and passed==total and total>0 and '[FAIL]' not in log and not diagnostic
            results.append(dict(program=wl.name,exit=code,passed=passed,total=total,
                                seconds=seconds,status='passed' if ok else 'failed'))
    receipt = dict(python=sys.version.split()[0],mathematica='requested' if args.mathematica else 'not run',results=results)
    (args.output/'results.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,indent=2))
    return 0 if all(r['status']=='passed' for r in results) else 1


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError,FileNotFoundError,subprocess.TimeoutExpired) as error:
        print('Verification failed:',error,file=sys.stderr)
        sys.exit(1)
