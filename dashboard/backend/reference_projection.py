"""Reference-only supplementary projection; no model rows or fitted axes change."""
import hashlib,json
import numpy as np
from scipy.spatial.distance import cdist

VERSION='revealed-reference-supplement-v1'
def sha_array(a):return hashlib.sha256(np.asarray(a,dtype='<f8').tobytes()).hexdigest()
def population_hash(ids):return hashlib.sha256(json.dumps(list(ids),separators=(',',':')).encode()).hexdigest()

def load_reference(path, expected_sha256, ligand):
    """Load a packet only through an explicit immutable file-hash binding."""
    from pathlib import Path
    p=Path(path)
    if any(q.is_symlink() for q in [p,*p.parents]):raise ValueError('Reference symlink forbidden')
    raw=p.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected_sha256:raise ValueError('Reference packet hash mismatch')
    ref=json.loads(raw)
    if ref['ligand_id']!=ligand:raise ValueError('Reference packet ligand mismatch')
    return ref

def cross_distances(reference, population, smiles, atoms=None):
    from backend.science import chemistry
    _,perms,_,_=chemistry(smiles)
    x=np.asarray(population,dtype=float);q=np.asarray(reference,dtype=float)
    if x.ndim!=3 or q.shape!=x.shape[1:] or not np.isfinite(x).all() or not np.isfinite(q).all():raise ValueError('Invalid reference/population coordinates')
    n=x.shape[1]
    if atoms is not None and any(type(a) is not int for a in atoms):raise ValueError('Atom selection must contain integer source indices')
    sel=np.arange(n) if atoms is None else np.array(sorted(set(atoms)),dtype=int)
    if not len(sel) or sel.min()<0 or sel.max()>=n:raise ValueError('Atom selection outside source graph')
    out=np.full(len(x),np.inf)
    for p in perms:
        if set(p[sel])!=set(sel):continue
        # The native symmetric kernel permutes the population (second argument).
        d=cdist(q[sel].reshape(1,-1),x[:,p[sel]].reshape(len(x),-1),'sqeuclidean')[0]/len(sel)
        np.minimum(out,d,out=out)
    return np.sqrt(out)

def supplementary(distance_matrix, saved_xy, cross_vector):
    """Gower cross-Gram contraction in the exact saved axes, without fitting.

    For Y=U sqrt(lambda), q=b U / sqrt(lambda). Y(Y'Y)^+ is the
    equivalent saved-coordinate operator, with zero axes retained as zero.
    """
    d=np.asarray(distance_matrix,dtype=float);y=np.asarray(saved_xy,dtype=float);v=np.asarray(cross_vector,dtype=float)
    n=len(d)
    if d.shape!=(n,n) or y.shape!=(n,2) or v.shape!=(n,) or not all(np.isfinite(a).all() for a in (d,y,v)):raise ValueError('Supplement shape/nonfinite mismatch')
    if np.any(v<0):raise ValueError('Negative reference distances')
    if not n:raise ValueError('Empty projection population')
    d2=d*d;v2=v*v
    b=-.5*(v2-d2.mean(axis=1)-v2.mean()+d2.mean())
    return b@y@np.linalg.pinv(y.T@y)

def reference_overlay(reference, analysis, distance_matrix, rows, coordinates, fe_distances, smiles):
    """Call after Store.analysis with its exact retained rows and distance array.

    Returns only a reference overlay object. Never alters analysis or its rows.
    The caller must pass rows after subset and metric-availability filtering.
    """
    spec=analysis['analysis_spec'];ids=analysis['pose_ids'];state=spec['state'];metric=spec['metric']
    if reference['ligand_id']!=spec['ligand'] or reference['analysis_snapshot']!=spec['snapshot']:raise ValueError('Reference ligand/snapshot mismatch')
    from rdkit import Chem
    if reference['graph'].get('status')=='AVAILABLE' and Chem.MolToSmiles(Chem.MolFromSmiles(smiles))!=reference['graph']['source_isomeric_smiles']:raise ValueError('Reference source graph mismatch')
    got=[r['id']+'::'+r['source_state'] if state=='paired' else r['id'] for r in rows]
    if got!=ids or population_hash(ids)!=spec['population_sha256']:raise ValueError('Reference population/order mismatch')
    d=np.asarray(distance_matrix,dtype=float)
    if sha_array(d)!=analysis['projection']['input_distance_sha256']:raise ValueError('Reference distance matrix binding mismatch')
    if sha_array(coordinates)!=analysis['coordinate_sha256']:raise ValueError('Reference coordinate binding mismatch')
    common=dict(id=reference['id'],ligand_id=reference['ligand_id'],role='revealed_reference',primary_copy=reference['primary_copy'],source_hashes=reference['source_hashes'],graph=reference['graph'],frame=reference['frame'],population_sha256=population_hash(ids),distance_sha256=sha_array(d),saved_xy_sha256=sha_array(analysis['xy']),coordinate_state=state,metric=metric,atom_indices=spec['atom_indices'],version=VERSION,refit=False)
    if metric=='contacts':return dict(common,status='UNAVAILABLE',xy=None,reason='Crystal contact fingerprint is not qualified; no RMSD coordinate reuse')
    if reference.get('coordinates') is None:return dict(common,status='UNAVAILABLE',xy=None,reason=reference['graph'].get('reason','Graph mapping unavailable'))
    if metric in ('rmsd','selected_atoms'):
        if metric=='selected_atoms' and not spec['atom_indices']:raise ValueError('Selected atoms required')
        if metric=='selected_atoms':
            choices=reference['graph'].get('alternative_source_to_native',[])
            selected=spec['atom_indices']
            if any(a<0 or a>=len(reference['coordinates']) for a in selected):raise ValueError('Selected atom outside source graph')
            sets={tuple(sorted(p[a] for a in selected)) for p in choices}
            if len(sets)>1:return dict(common,status='UNAVAILABLE',xy=None,reason='Selected source atom set is not invariant across crystal graph symmetries; no unique crystal atom correspondence')
        v=cross_distances(reference['coordinates'],coordinates,smiles,spec['atom_indices'] if metric=='selected_atoms' else None)
    elif metric=='heme':
        if reference.get('fe_distances') is None:return dict(common,status='UNAVAILABLE',xy=None,reason='Primary heme Fe unavailable')
        v=cross_distances(np.asarray(reference['fe_distances'])[:,None],np.asarray(fe_distances)[:,:,None],smiles)
    else:return dict(common,status='UNAVAILABLE',xy=None,reason='Unsupported reference metric')
    xy=supplementary(d,analysis['xy'],v)
    return dict(common,status='AVAILABLE',xy=xy.tolist(),cross_distance_sha256=sha_array(v),reason=None,method='Gower supplementary PCoA in saved axes; no full n+1 matrix')
