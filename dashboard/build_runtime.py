#!/usr/bin/env python3
"""Build an isolated release from this checkout and operator-provided data."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / 'runtime-data.manifest.json'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def safe_path(root, relative):
    name = PurePosixPath(relative)
    if name.is_absolute() or '..' in name.parts or not name.parts:
        raise ValueError('Invalid relative path: ' + relative)
    path = root / relative
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError('Path leaves supplied root: ' + relative)
    return path


def verify_file(path, record, hash_key='sha256', size_key='bytes'):
    if not path.is_file() or path.stat().st_size != record[size_key] or digest(path) != record[hash_key]:
        raise ValueError('Missing or changed input: ' + str(path))


def verify_coordinates(root, coordinate_manifest, manifest):
    if digest(coordinate_manifest) != manifest['coordinate_manifest_sha256']:
        raise ValueError('Coordinate manifest hash mismatch')
    records = json.loads(coordinate_manifest.read_text())['files']
    if len(records) != manifest['coordinate_count'] or len({r['path'] for r in records}) != len(records):
        raise ValueError('Coordinate membership mismatch')
    for row in records:
        verify_file(safe_path(root, row['path']), row, 'gzip_sha256', 'gzip_bytes')
    return records


def build(data_root, coordinate_root, coordinate_manifest, output):
    manifest = json.loads(MANIFEST.read_text())
    output = output.absolute()
    if output.exists() or output.is_symlink():
        raise ValueError('Output must not exist')
    if output.resolve().is_relative_to(ROOT):
        raise ValueError('Build outside the canonical source directory')
    data_root = data_root.resolve()
    coordinate_root = coordinate_root.resolve()
    records = verify_coordinates(coordinate_root, coordinate_manifest, manifest)
    for name, row in manifest['files'].items():
        verify_file(safe_path(data_root, name), row)
    for name, row in manifest['code_files'].items():
        verify_file(safe_path(ROOT, name), row)
    subprocess.run([sys.executable, str(ROOT / 'build_assets.py'), '--check'], check=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.' + output.name + '-', dir=output.parent))
    try:
        for collection, source in [(manifest['files'], data_root), (manifest['code_files'], ROOT)]:
            for name, row in collection.items():
                target = safe_path(stage, name)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(safe_path(source, name), target)
                verify_file(target, row)
        for name in ['publication_entry.py', 'run_runtime.py', 'verify_runtime.py']:
            shutil.copyfile(ROOT / name, stage / name)
        shutil.copytree(ROOT / 'licenses', stage / 'licenses')
        # Bind external coordinates to the structure allowlist, not just a count.
        structures = json.loads((stage / 'structures.json').read_text())
        lookup = {r['path']: r for r in records}
        if {v['path'] for v in structures.values()} != set(lookup):
            raise ValueError('Structure allowlist and coordinate manifest differ')
        for value in structures.values():
            if value['sha256'] != lookup[value['path']]['sha256']:
                raise ValueError('Uncompressed coordinate hash mismatch')
        env = dict(os.environ, CYP_RECEPTOR31_DIR=str(stage / 'data/receptor-map31'),
                   CYP_COORDINATE_ROOT=str(coordinate_root), PYTHONPATH=str(stage),
                   PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        probe = subprocess.check_output([sys.executable, '-c', 'from backend.runtime import SID; print(SID)'], cwd=stage, env=env, text=True).strip()
        if probe != manifest['runtime_snapshot']:
            raise ValueError('Reconstructed runtime snapshot mismatch: ' + probe)
        try:
            commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
            dirty = bool(subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=normal'], cwd=ROOT, text=True))
        except subprocess.CalledProcessError:
            commit, dirty = None, True
        receipt = {'schema': 1, 'source_commit': commit, 'source_dirty': dirty,
                   'runtime_snapshot': probe, 'runtime_data_manifest_sha256': digest(MANIFEST),
                   'runtime_files': len(manifest['files']) + len(manifest['code_files']),
                   'coordinate_files_verified': len(records),
                   'files': {str(p.relative_to(stage)): {'sha256': digest(p), 'bytes': p.stat().st_size}
                             for p in sorted(stage.rglob('*')) if p.is_file() and 'cache' not in p.relative_to(stage).parts}}
        (stage / 'build-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        stage.rename(output)
        print(json.dumps({k: v for k, v in receipt.items() if k != 'files'}, indent=2))
    except BaseException:
        shutil.rmtree(stage)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', required=True, type=Path)
    parser.add_argument('--coordinate-root', required=True, type=Path)
    parser.add_argument('--coordinate-manifest', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    build(args.data_root, args.coordinate_root, args.coordinate_manifest, args.output)


if __name__ == '__main__':
    main()
