"""Presentation universe only; immutable scientific inputs remain untouched."""
import json
from collections import Counter

VERSION = 'dashboard-universe-v3-ipsae2-water'
EXCLUDED = frozenset({'boltz-metadiffusion-off', 'flowdock'})

def included(method):
    return method not in EXCLUDED

def configure(store, data):
    store.ledger = {k:v for k,v in store.ledger.items() if included(v['method'])}
    store.receptor = {k:v for k,v in store.receptor.items() if k.split('|',1)[-1] in store.ledger}
    store.normalization = dict(store.normalization)
    store.normalization['groups'] = [g for g in store.normalization.get('groups',[]) if included(g['method'])]
    counts_path=data/'dashboard_counts.json'
    frozen_counts=json.loads(counts_path.read_text()) if counts_path.exists() else None
    raw = relaxed = 0
    for lid, info in store.ligands.items():
        counts = {}
        for state in ('raw','relaxed'):
            if frozen_counts is not None:counts[state]=frozen_counts[lid][state]
            else:
                from .historical import feature_folder
                rows = json.loads((feature_folder(data,lid,state)/(state+'.json')).read_text())
                counts[state] = sum(included(r['method']) for r in rows)
        info.update(n_total=sum(r['ligand_id']==lid for r in store.ledger.values()), n_raw=counts['raw'], n_relaxed=counts['relaxed'])
        raw += counts['raw']; relaxed += counts['relaxed']
    store.snapshot = dict(store.snapshot, id=store.snapshot['id']+':'+VERSION,
        total=len(store.ledger), verified_relaxed=relaxed,
        statuses=dict(Counter(r['status'] for r in store.ledger.values())),
        scope='Frozen source snapshot; dashboard presentation universe')
    store.feature_coverage = dict(ligands=len(store.ligands), raw=raw, relaxed=relaxed,
        receptor_rows=len(store.receptor), unavailable=len(store.ledger)-raw,
        scope='Dashboard presentation universe; no new refinement or inference')
