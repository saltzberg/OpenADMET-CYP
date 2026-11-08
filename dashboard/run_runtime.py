#!/usr/bin/env python3
"""Start one CPU-only worker from a verified release."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coordinate-root', required=True, type=Path)
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', type=int, default=8896)
    args = parser.parse_args()
    receipt = json.loads((ROOT / 'build-receipt.json').read_text())
    for name, record in receipt['files'].items():
        target = ROOT / name
        if not target.resolve().is_relative_to(ROOT) or not target.is_file():
            raise ValueError('Release file is missing or outside the release: ' + name)
        payload = target.read_bytes()
        if len(payload) != record['bytes'] or hashlib.sha256(payload).hexdigest() != record['sha256']:
            raise ValueError('Release file changed: ' + name)
    if receipt['source_commit'] and not receipt['source_dirty']:
        os.environ['CYP_SOURCE_COMMIT'] = receipt['source_commit']
    else:
        os.environ.pop('CYP_SOURCE_COMMIT', None)
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[name] = '1'
    os.environ.update(CYP_RECEPTOR31_DIR=str((ROOT / 'data/receptor-map31').resolve()),
                      CYP_COORDINATE_ROOT=str(args.coordinate_root.resolve()),
                      CUDA_VISIBLE_DEVICES='', PYTHONDONTWRITEBYTECODE='1', PYTHONPATH=str(ROOT))
    os.chdir(ROOT)
    os.execv(sys.executable, [sys.executable, '-m', 'uvicorn', 'publication_entry:app',
                            '--host', args.host, '--port', str(args.port), '--workers', '1',
                            '--limit-concurrency', '24', '--backlog', '24',
                            '--timeout-keep-alive', '3', '--timeout-graceful-shutdown', '8'])


if __name__ == '__main__':
    main()
