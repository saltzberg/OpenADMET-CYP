"""Optional frozen historical coordinate overlay. No scientific computation here."""
import json
import pandas as pd

def configure(store):
    store.historical_energy = {}
    for row in store.read('enrichment_historical_ledger.json', []):
        pid = row['pose_id']
        assert pid in store.ledger and store.ledger[pid]['raw_sha256'] == row['raw_sha256']
        store.ledger[pid] = row
    store.receptor.update(store.read('enrichment_historical_receptor.json', {}))
    path = store.data/'enrichment_historical_energies.parquet'
    if path.exists():
        store.historical_energy = pd.read_parquet(path).set_index('pose_id').to_dict('index')
    store.historical = store.read('enrichment_historical_build.json', {})
    if store.historical:
        store.snapshot = dict(store.snapshot, id=store.snapshot['id']+':historical-'+store.historical['source_sha256'][:12], historical_relaxed=store.historical['notice'], historical_selected=store.historical['selected_count'], historical_added=store.historical['new_count'])

def feature_folder(data, ligand, state):
    folder = data/'historical_features'/ligand
    return folder if state == 'relaxed' and (folder/'relaxed.json').exists() else data/'features'/ligand

def energies(store, values):
    for pid, row in store.historical_energy.items():
        source = store.ledger[pid]
        assert row['selected_result_sha256'] == source['result_sha256']
        assert row['selected_structure_sha256'] == source['source_hashes'][source['relaxed_structure_path']]
        assert row['selected_positions_sha256'] == source['source_hashes'][source['coordinates_path']]
        values[pid] = dict(row, pose_id=pid)
    return values
