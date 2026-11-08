"""Source-graph geometry and conservative typed contacts (not ProLIF)."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.spatial.distance import cdist,pdist,squareform
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.linalg import eigh
from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures
from functools import lru_cache
import json
import hashlib

PROJECTION_VERSION='pcoa-canonical-eigenspaces-v2'
CLUSTER_ALGORITHM='average-linkage-original-distances-v1'

def source_graph_coordinates(x, native, source_smiles, source_to_native):
 """Validate explicit full graph/stereo correspondence; never fit coordinates.

 Native must carry authoritative bond orders/charge/stereo, not distance-inferred
 bonds. Unspecified source stereo is not upgraded to coordinate stereo evidence.
 """
 source=Chem.MolFromSmiles(source_smiles)
 if source is None or native is None:raise ValueError('invalid ligand graph')
 n=source.GetNumAtoms();p=tuple(source_to_native)
 if native.GetNumAtoms()!=n or native.GetNumBonds()!=source.GetNumBonds() or sorted(p)!=list(range(n)):
  raise ValueError('not a full graph bijection')
 matches=native.GetSubstructMatches(source,uniquify=False,useChirality=True,maxMatches=100001)
 if len(matches)>100000:raise ValueError('graph mapping enumeration exceeds bound')
 # Bidirectional isomeric identity rejects state/stereo specified on only one side.
 if Chem.MolToSmiles(native,isomericSmiles=True)!=Chem.MolToSmiles(source,isomericSmiles=True) or p not in matches:
  raise ValueError('graph/state/stereo correspondence mismatch')
 x=np.asarray(x,dtype=float)
 if x.shape!=(n,3) or not np.isfinite(x).all():raise ValueError('invalid native coordinates')
 return x[list(p)].copy()

def project_receptor(features, mean, basis):
 """Transform only through a caller-supplied frozen, versioned PCA basis."""
 x=np.asarray(features,dtype=float);mean=np.asarray(mean,dtype=float);basis=np.asarray(basis,dtype=float)
 if x.ndim!=2 or mean.shape!=(x.shape[1],) or basis.shape!=(x.shape[1],2):raise ValueError('receptor basis shape mismatch')
 if not all(np.isfinite(a).all() for a in (x,mean,basis)):raise ValueError('nonfinite receptor features/basis')
 # Rowwise contraction makes existing points independent of batch length.
 return np.einsum('ij,jk->ik',x-mean,basis,optimize=False)

def _canonical_pcoa(d):
 """Top positive eigenspaces, with index-anchored degenerate orientations.

 Expand the solved range if a tie crosses its boundary. Anchor a repeated
 eigenspace via its projector onto coordinate unit vectors in input row order;
 this is invariant to the arbitrary eigenvectors returned by LAPACK. Numerical
 ties use 1e-10 relative spectral tolerance. Input ordering is part of identity.
 """
 n=len(d)
 if not np.any(d):return np.zeros((n,2))
 b=-.5*d*d;b-=b.mean(0)[None,:];b-=b.mean(1)[:,None]
 k=min(n,8)
 while True:
  v,u=eigh(b,subset_by_index=[n-k,n-1],driver='evr');v=v[::-1];u=u[:,::-1]
  tol=max(np.finfo(float).tiny,float(abs(v).max()))*1e-10
  if k==n or abs(v[1]-v[-1])>tol:break
  k=min(n,k*2)
 columns=[];start=0
 while len(columns)<2 and start<len(v) and v[start]>tol:
  stop=start+1
  while stop<len(v) and abs(v[stop]-v[start])<=tol:stop+=1
  block=u[:,start:stop]
  if stop-start==1:
   q=block[:,0].copy()
   if q[np.argmax(abs(q))]<0:q*=-1
   columns.append(q*np.sqrt(v[start]))
  else:
   anchors=[]
   for i in range(n):
    q=block@block[i,:]
    for prev in anchors:q-=prev*np.dot(prev,q)
    norm=np.linalg.norm(q)
    if norm>1e-9:
     q/=norm;anchors.append(q);columns.append(q*np.sqrt(float(v[start:stop].mean())))
     if len(columns)==2:break
  start=stop
 return np.column_stack(columns+[np.zeros(n)]*(2-len(columns)))
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
SETS=json.loads((BASE/'pocket-residue-sets-20261028/residue_sets.json').read_text())
@lru_cache(maxsize=64)
def chemistry(smiles):
 m=Chem.MolFromSmiles(smiles)
 if m is None:raise ValueError('invalid source SMILES')
 perms=m.GetSubstructMatches(m,uniquify=False,useChirality=True,maxMatches=100001)
 if len(perms)>100000:raise ValueError('symmetry enumeration exceeds bound')
 # Full graph matches preserve stereo; no element-only assignment.
 ranks=list(Chem.CanonicalRankAtoms(m,breakTies=False,includeChirality=True))
 factory=ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef'))
 features=factory.GetFeaturesForMol(m)
 types={k:set() for k in ('Donor','Acceptor','Hydrophobe','PosIonizable','NegIonizable')}
 for f in features:
  if f.GetFamily() in types:types[f.GetFamily()].update(f.GetAtomIds())
 return m,np.array(perms,dtype=int),ranks,types

def rmsd_matrix(x,smiles,atoms=None):
 m,perms,_,_=chemistry(smiles); n=m.GetNumAtoms()
 x=np.asarray(x,dtype=float)
 if len(x) and (x.ndim!=3 or x.shape[1]!=n or x.shape[2] not in (1,3) or not np.isfinite(x).all()):raise ValueError('expected finite source-ordered coordinates (poses, atoms, 3) or Fe distances (poses, atoms, 1)')
 if atoms is not None and any(not isinstance(a,(int,np.integer)) or isinstance(a,(bool,np.bool_)) for a in atoms):raise ValueError('atom selection must be integer source indices')
 selected=np.arange(n) if atoms is None else np.array(sorted(set(atoms)),dtype=int)
 if len(selected)==0 or selected.min()<0 or selected.max()>=n:raise ValueError('atom selection outside source graph')
 if len(x)==0:return np.zeros((0,0))
 perms=np.array([p for p in perms if set(p[selected])==set(selected)])
 d=np.full((len(x),len(x)),np.inf)
 a=x[:,selected].reshape(len(x),-1)
 for p in perms:
  q=cdist(a,x[:,p[selected]].reshape(len(x),-1),'sqeuclidean')/len(selected)
  np.minimum(d,q,out=d)
 np.sqrt(d,out=d);d=(d+d.T)/2;np.fill_diagonal(d,0)
 return d

def contact_matrix(rows):
 keys=sorted(set(k for r in rows for k in r));look={k:i for i,k in enumerate(keys)}
 from scipy.sparse import csr_matrix
 ii=[];jj=[]
 for i,r in enumerate(rows):
  for k in set(r):ii.append(i);jj.append(look[k])
 a=csr_matrix((np.ones(len(ii)),(ii,jj)),shape=(len(rows),len(keys)))
 dot=(a@a.T).toarray();counts=np.asarray(a.sum(1)).ravel();union=counts[:,None]+counts[None,:]-dot
 d=1-np.divide(dot,union,out=np.zeros_like(dot),where=union>0)
 # Empty signatures are not evidence of a shared bound mode.
 np.fill_diagonal(d,0);return d

def analyze_distances(d,ids,threshold,metric,*,precomputed=None):
 n=len(ids);d=np.asarray(d,dtype=float)
 if d.shape!=(n,n) or not np.isfinite(d).all() or np.any(d<0) or not np.allclose(d,d.T,rtol=0,atol=1e-10) or np.any(np.diag(d)!=0):raise ValueError('invalid distance matrix')
 if len(set(ids))!=n:raise ValueError('duplicate projection pose IDs')
 if not np.isfinite(threshold) or threshold<0:raise ValueError('threshold must be finite and nonnegative')
 if precomputed is not None:
  xy=np.asarray(precomputed['xy'],dtype=float);tree=np.asarray(precomputed['linkage'],dtype=float)
  if xy.shape!=(n,2) or tree.shape!=(max(0,n-1),4) or not np.isfinite(xy).all() or not np.isfinite(tree).all():raise ValueError('invalid precomputed analysis state')
 projection=dict(method='Deterministic classical metric MDS / PCoA',metric=metric,version=PROJECTION_VERSION,deterministic=True,seed=None,orientation='input-row anchored eigenspaces; numerical tie tolerance 1e-10 relative; same ordered inputs and numerical runtime',distortion_label='2D normalized distance stress',distortion_definition='sqrt(sum((D - d2D)^2) / sum(D^2)); no fitted rescaling',distortion_status='available',cluster_space='original distances; average linkage',input_distance_sha256=hashlib.sha256(np.asarray(d,dtype='<f8').tobytes()).hexdigest())
 if n==0:
  projection.update(stress=None,distortion_status='unavailable_no_points')
  return dict(pose_ids=[],xy=[],labels=[],clusters=[],projection=projection,warnings=['No coordinates available'])
 tree=precomputed['linkage'] if precomputed is not None else linkage(squareform(d,checks=False),method='average') if n>1 else None
 labels=fcluster(tree,threshold,criterion='distance') if n>1 else np.ones(1,int)
 if precomputed is not None:xy=np.asarray(precomputed['xy'],dtype=float)
 elif n>2:xy=_canonical_pcoa(d)
 elif n==2:xy=np.array([[-d[0,1]/2,0],[d[0,1]/2,0]])
 else:xy=np.zeros((1,2))
 if precomputed is not None and 'projection' in precomputed:
  projection=dict(precomputed['projection'])
 else:
  den=float(np.sum(d*d));stress=float(np.sqrt(np.sum((d-squareform(pdist(xy)))**2)/den)) if den else 0.
  projection.update(stress=stress,distortion_status='available' if den else 'zero_distance_convention')
 groups=[]
 for label in sorted(set(labels)):
  ix=np.flatnonzero(labels==label);sub=d[np.ix_(ix,ix)];rep=int(ix[np.argmin(sub.mean(1))]);chosen=[rep]
  while len(chosen)<min(5,len(ix)):
   remaining=[int(i) for i in ix if i not in chosen];chosen.append(max(remaining,key=lambda i:(float(d[i,chosen].min()),-i)))
  groups.append(dict(id=int(label),count=len(ix),representative_id=ids[rep],member_ids=[ids[i] for i in ix],display_ids=[ids[i] for i in chosen]))
 return dict(pose_ids=ids,xy=xy.tolist(),labels=labels.tolist(),clusters=groups,projection=projection,warnings=['2D projection is diagnostic, not clustering space.'])

# Exact canonical identities; equivalence only for genuine residue symmetries.
PROTEIN_EQ={'ASP':{'OD1':'OD*','OD2':'OD*'},'GLU':{'OE1':'OE*','OE2':'OE*'},'ARG':{'NH1':'NH*','NH2':'NH*'},'PHE':{'CD1':'CD*','CD2':'CD*','CE1':'CE*','CE2':'CE*'},'TYR':{'CD1':'CD*','CD2':'CD*','CE1':'CE*','CE2':'CE*'},'VAL':{'CG1':'CG*','CG2':'CG*'},'LEU':{'CD1':'CD*','CD2':'CD*'}}
def typed_contacts(x,env,smiles,cutoff=4.5):
 m,perms,ranks,types=chemistry(smiles);xyz=np.array([a['xyz'] for a in env]);dist=cdist(x,xyz);records=[];keys=set()
 for i,j in zip(*np.where(dist<=cutoff)):
  a=env[j];e=a['element'].upper();res=a['resname'];name=a['name'];d=float(dist[i,j]);kind=None
  if a['chain']=='H':
   if e=='FE':kind='Fe_donor_proximity' if i in types['Acceptor'] and d<=3.0 else 'Fe_proximity'
   else:kind='heme_atom_proximity'
  elif a['chain']=='A':
   if m.GetAtomWithIdx(int(i)).GetSymbol()=='C' and e=='C' and i in types['Hydrophobe'] and name not in ('C','CA'):kind='hydrophobic_proximity'
   # Polar candidates deliberately do not assert H bonds without H-angle proof.
   elif ((i in types['Donor'] and e in ('O','N','S')) or (i in types['Acceptor'] and e in ('N','O','S'))) and d<=3.5:kind='polar_contact_candidate'
   else:kind='heavy_atom_proximity'
  if kind is None:continue
  protein_name=PROTEIN_EQ.get(res,{}).get(name,name)
  site=f'{a["chain"]}:{res}{a["resid"]}:{protein_name}'
  if a['chain']=='H' and e=='FE':site='H:HEM:FE'
  key=f'Lsym{ranks[i]}|{site}|{kind}'
  # Non-Fe native heme atom names differ across methods. Keep exact atom detail,
  # but do not manufacture a cross-method correspondence in the fingerprint.
  fingerprint_member=not (a['chain']=='H' and e!='FE')
  if fingerprint_member:keys.add(key)
  records.append(dict(key=key,fingerprint_member=fingerprint_member,ligand_atom=int(i),ligand_symmetry_class=int(ranks[i]),protein_chain=a['chain'],protein_residue=res,protein_position=a['resid'],protein_atom=name,interaction_type=kind,distance_A=d))
 return sorted(keys),records

def receptor_features(env,residues):
 ca={int(a['resid']):np.array(a['xyz']) for a in env if a['chain']=='A' and a['name']=='CA' and a['resid'] is not None}
 if not all(r in ca for r in residues):return None
 centers=[]
 for r in residues:
  side=[a['xyz'] for a in env if a['chain']=='A' and a['resid'] is not None and int(a['resid'])==r and a['name'] not in ('N','CA','C','O','OXT') and a['element'].upper()!='H']
  centers.append(np.mean(side,axis=0) if side else ca[r])
 anchors=np.array([ca[r] for r in (residues[0],residues[len(residues)//3],residues[2*len(residues)//3],residues[-1])])
 return np.concatenate([cdist(np.array([ca[r] for r in residues]),anchors).ravel(),cdist(np.array(centers),anchors).ravel()])
