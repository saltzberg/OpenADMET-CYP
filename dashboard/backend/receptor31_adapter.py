"""Opt-in adapter. This module does not edit packets or serving files.

Call configure(store, artifact_folder) after FrozenStore construction.
Include scientific_revision(folder) in the outer runtime snapshot/cache ID.
Keep store.snapshot unchanged: it identifies frozen pose/energy observations.
"""
import copy,hashlib,json,types
from pathlib import Path
import numpy as np

def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def _load(folder):
    folder=Path(folder)
    hashes=json.loads((folder/'artifact-hashes.json').read_text())
    for name in ('manifest.json','receptor_basis.npz','receptor_map.json','revealed_map.json'):
        if _sha(folder/name)!=hashes[name]:raise ValueError('Receptor31 artifact hash mismatch: '+name)
    manifest=json.loads((folder/'manifest.json').read_text())
    poses=json.loads((folder/'receptor_map.json').read_text())
    refs=json.loads((folder/'revealed_map.json').read_text())
    if manifest['basis_sha256']!=_sha(folder/'receptor_basis.npz'):raise ValueError('Basis identity mismatch')
    if len(poses)!=manifest['count'] or len(refs)!=20:raise ValueError('Coverage mismatch')
    with np.load(folder/'receptor_basis.npz',allow_pickle=False) as z:
        if z['mean'].shape!=(248,) or z['basis'].shape!=(248,2) or list(z['canonical_residues'])!=manifest['canonical_residues']:raise ValueError('Basis schema mismatch')
    if not np.isfinite(np.array(list(poses.values()))).all():raise ValueError('Nonfinite receptor map')
    return manifest,poses,refs

def scientific_revision(folder):
    """Bind response caches, background caches and runtime fingerprint to map."""
    folder=Path(folder)
    _load(folder)
    return hashlib.sha256(b'receptor31-adapter-v1\0'+b''.join((folder/n).read_bytes() for n in ('manifest.json','receptor_basis.npz','receptor_map.json','revealed_map.json'))).hexdigest()

def decorate_bundle(bundle,refs,revision):
    """Replace only receptor overlay. Preserve pose maps, energies and selectors."""
    out=copy.deepcopy(bundle)
    ref=refs[bundle['ligand']]
    if ref['source_sha256']!=bundle['source_pdb_sha256'] or ref['primary_copy']!=bundle['primary_copy'] or ref['id']!=bundle['id']:raise ValueError('Revealed source/copy mismatch')
    out['receptor']={**ref,'method':'Shared 31-residue receptor PCA; supplementary transform only','refit':False,'reason':None}
    out['receptor_map_revision']=revision
    out['revision']=hashlib.sha256((str(bundle.get('revision',''))+'|'+revision).encode()).hexdigest()
    return out

def configure(store,folder):
    """Install an exact-ID receptor map; all other scientific data stay intact."""
    if getattr(store,'_receptor31_configured',False):raise ValueError('Receptor31 already configured')
    manifest,poses,refs=_load(folder)
    if set(poses)!=set(store.receptor):raise ValueError('Receptor31 does not match current frozen pose universe')
    revision=scientific_revision(folder)
    store.receptor=poses
    store.projection_version=manifest['projection_version']+':'+revision
    store.receptor_projection={**manifest,'revision':revision,'scientific_source_snapshot':getattr(store,'source_snapshot',None)}
    original_analysis=store.analysis
    def analysis(self,*args,**kwargs):
        result=original_analysis(*args,**kwargs)
        if 'revealed_reference' in result:
            result=dict(result,revealed_reference=decorate_bundle(result['revealed_reference'],refs,revision))
            self._revealed_reference=result['revealed_reference']
        result=dict(result,receptor_projection=self.receptor_projection,receptor_map_revision=revision)
        return result
    store.analysis=types.MethodType(analysis,store)
    store._receptor31_configured=True
    return revision
