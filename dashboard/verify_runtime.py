#!/usr/bin/env python3
"""Check release bytes and real HTTP data. Requires provisioned scientific inputs."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlopen

PREFIX = '/articles/admet-cyp-structure/cofolding/'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime-root', type=Path, required=True)
    parser.add_argument('--origin', required=True, help='Origin only, for example http://127.0.0.1:8896')
    parser.add_argument('--all-ligands', action='store_true')
    parser.add_argument('--source-commit', help='Require this X-Source-Commit header')
    args = parser.parse_args()
    root = args.runtime_root.resolve()
    sys.path.insert(0, str(root))
    from backend.compact import decode
    receipt = json.loads((root / 'build-receipt.json').read_text())
    for name, row in receipt['files'].items():
        data = (root / name).read_bytes()
        assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'], name
    def get(path):
        with urlopen(args.origin.rstrip('/') + path, timeout=120) as response:
            if args.source_commit:
                assert response.headers.get('X-Source-Commit') == args.source_commit
            return response.read()
    boot = json.loads(get(PREFIX + 'api/bootstrap'))
    assert boot['snapshot'] == receipt['runtime_snapshot']
    html = get(PREFIX).decode()
    assert '?v=' + boot['snapshot'] in html
    base = boot['api_base']
    meta = json.loads(get(base + '/meta'))
    projection = json.loads((root / 'data/receptor-map31/manifest.json').read_text())
    assert meta['receptor_projection']['projection_version'] == projection['projection_version']
    receptor = json.loads((root / 'data/receptor-map31/receptor_map.json').read_text())
    structures = json.loads((root / 'structures.json').read_text())
    lids = [row['id'] for row in meta['ligands']] if args.all_ligands else [meta['default_ligand']]
    counts = {'raw': 0, 'relaxed': 0, 'paired': 0}
    scored = {'raw': 0, 'relaxed': 0}
    for lid in lids:
        source = json.loads((root / 'data/truth-metrics' / (lid + '.json')).read_text())
        truth = {(r['source_state'], r['pose_id']): r for r in source}
        for state in ('raw', 'relaxed', 'paired'):
            wire = json.loads(get(base + '/ligand?ligand=' + lid + '&state=' + state + '&format=columns-v1'))
            assert 'columns' in wire and wire['snapshot'] == boot['snapshot']
            poses = decode(wire)['poses']
            counts[state] += len(poses)
            for pose in poses:
                native = pose['id'].removesuffix('::' + pose['source_state'])
                key = (pose['source_state'], native)
                for metric in ('LDDT-PLI', 'BiSyRMSD'):
                    assert pose['truth_metrics'][metric] == truth[key][metric], (lid, key, metric)
                assert pose['structure_sha256'] == truth[key]['structure_sha256']
                assert pose['receptor_xy'] == receptor['|'.join(key)], (lid, key, 'map')
                if state != 'paired' and pose['truth_metrics']['LDDT-PLI'] is not None:
                    scored[state] += 1
            first = next(p for p in poses if p.get('structure_sha256'))
            match = next(v for v in structures.values() if v['sha256'] == first['structure_sha256'])
            assert hashlib.sha256(gzip.decompress(get(PREFIX + match['path']))).hexdigest() == match['sha256']
            print(json.dumps({'ligand': lid, 'state': state, 'rows': len(poses), 'metrics_and_map': 'PASS', 'coordinate': 'PASS'}), flush=True)
    if args.all_ligands:
        assert counts == {'raw': 22000, 'relaxed': 20213, 'paired': 40426}, counts
        assert scored == {'raw': 21934, 'relaxed': 20213}, scored
    print(json.dumps({'status': 'PASS', 'snapshot': boot['snapshot'], 'ligands': len(lids), 'counts': counts, 'scored': scored}))


if __name__ == '__main__':
    main()
