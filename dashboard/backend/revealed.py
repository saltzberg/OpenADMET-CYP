"""Hash-bound experimental reference overlays; never generated-model rows."""
import hashlib,json
from pathlib import Path
import numpy as np
from .reference_projection import load_reference,reference_overlay

VERSION='revealed-reference-v1'

def file_sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_packet(store,ligand,folder=None):
    folder=Path(folder) if folder else store.data/'revealed-reference'
    path=folder/'manifest.json'
    if not path.is_file():raise ValueError('Revealed reference manifest unavailable')
    manifest=json.loads(path.read_text())
    if manifest['version']!=VERSION:raise ValueError('Unknown reference version')
    source=json.loads((store.data.parent/'manifest.json').read_text())['snapshot']
    if manifest['source_snapshot']!=source or manifest['analysis_snapshot']!=store.snapshot['id']:raise ValueError('Reference source snapshot mismatch')
    item=manifest['ligands'][ligand]
    for name,digest in item['feature_sources'].items():
        target=store.data/name
        if '..' in Path(name).parts or Path(name).is_absolute():raise ValueError('Invalid reference feature path')
        if file_sha(target)!=digest:raise ValueError('Reference feature source mismatch')
    name=item['file']
    if Path(name).name!=name:raise ValueError('Invalid reference packet path')
    packet=load_reference(folder/name,item['sha256'],ligand)
    if packet['analysis_snapshot']!=manifest['analysis_snapshot']:raise ValueError('Reference analysis snapshot mismatch')
    if packet['source_pdb_sha256']!=item['source_pdb_sha256']:raise ValueError('Reference input hash mismatch')
    return packet,file_sha(path)

def decorate_analysis(store,result,folder=None):
    """Append only an overlay, bound to the exact saved ordered analysis."""
    spec=result['analysis_spec'];ligand=spec['ligand'];state=spec['state']
    reference,revision=load_packet(store,ligand,folder)

    rows,x,fe=store.features(ligand,state)
    wanted=set(result['pose_ids'])
    ix=[i for i,r in enumerate(rows) if (r['id']+'::'+r['source_state'] if state=='paired' else r['id']) in wanted]
    rows=[rows[i] for i in ix];x=x[ix];fe=fe[ix]
    key=result['diagnostics']['phase_cache_key']
    matrix=np.load(store.analysis_folder/(key+'.npy'),allow_pickle=False)
    if len(rows):pose=reference_overlay(reference,result,matrix,rows,x,fe,store.ligands[ligand]['smiles'])
    else:pose=dict(id=reference['id'],role='revealed_reference',status='UNAVAILABLE',xy=None,reason='No retained model coordinates for supplementary projection',metric=spec['metric'],coordinate_state=state,atom_indices=spec['atom_indices'],refit=False)
    energy=reference['energy']
    selected=[r for r in store.features(ligand,'relaxed')[0] if r.get('final_selection') and r['source_state']=='relaxed']
    if len(selected)!=1 or selected[0]['structure_sha256']!=energy['submitted_structure_sha256']:raise ValueError('Reference submitted structure mismatch')
    # Even an unavailable energy is bound to its intended exact row reference.
    old=store.energies().get(selected[0]['id'],{})
    for field in ('reference_key','reference_version','protocol_id'):
        if old.get(field)!=energy['submitted_'+field]:raise ValueError('Reference isolated-energy binding mismatch')
    bundle=dict(version=VERSION,revision=revision,snapshot=result['snapshot'],source_snapshot=store.source_snapshot,analysis_snapshot=reference['analysis_snapshot'],ligand=ligand,state=state,id=reference['id'],role='revealed_reference',primary_copy=reference['primary_copy'],secondary_copies=reference['secondary_copies'],source_pdb_sha256=reference['source_pdb_sha256'],source_hashes=reference['source_hashes'],graph=reference['graph'],observed_receptor=reference['observed_receptor'],pose=pose,receptor=reference['receptor_projection'],energy=energy)
    result=dict(result,revealed_reference=bundle)
    store._revealed_reference=bundle
    return result
