"""Scientific store and legacy in-process test adapter.
Production: .venv-service/bin/python -m backend.main (FastAPI/uvicorn only).
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,hashlib,argparse,mimetypes,threading,math
from pathlib import Path
from functools import lru_cache
from urllib.parse import urlparse,parse_qs
import numpy as np,pandas as pd,gemmi
from .science import rmsd_matrix,contact_matrix,analyze_distances,PROJECTION_VERSION,CLUSTER_ALGORITHM
from .dashboard_universe import VERSION, included, configure
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data';COMPUTE=threading.BoundedSemaphore(1)
METRIC_LABELS={'pl_iptm':'Native protein → ligand iPTM','pl_iptm_reverse':'Native ligand → protein iPTM','heme_ligand_iptm':'Native heme → ligand iPTM','ligand_heme_iptm':'Native ligand → heme iPTM','ligand_local_confidence_native':'Native ligand local confidence','global_local_confidence_native':'Native global local confidence','global_iptm':'Native global iPTM','ipsae__L_to_Pmedium__row_coverage_mean':'ipSAE-like ligand → medium pocket mean','ipsae__Pmedium_to_L__row_coverage_mean':'ipSAE-like medium pocket → ligand mean','ipsae__L_to_H__row_coverage_mean':'ipSAE-like ligand → heme mean'}
PALETTE=['#4c78a8','#72b7b2','#54a24b','#e45756','#f58518','#b279a2','#9d755d','#bab0ac','#d4a72c','#5f6b8c','#b95a84','#597b58','#707070']
def read(name,default=None):
 p=DATA/name
 return json.loads(p.read_text()) if p.exists() else default

def clean(x):
 if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v) for v in x]
 if isinstance(x,np.ndarray):return clean(x.tolist())
 if isinstance(x,np.generic):return clean(x.item())
 if isinstance(x,float) and not math.isfinite(x):return None
 return x

class Store:
 def __init__(self,data=None,source_manifest=None):
  self.data=Path(data) if data else DATA
  self.read=lambda name,default=None: json.loads((self.data/name).read_text()) if (self.data/name).exists() else default
  self.source_manifest=Path(source_manifest) if source_manifest else ROOT.parent/'restrained-minimization-20261027/source_manifest.jsonl'
  self.full_features=(self.data/'features_complete.json').exists()
  self.snapshot=self.read('snapshot.json',self.read('snapshot_preview.json',{}));self.ligands={r['id']:r for r in self.read('ligands.json' if self.full_features else 'ligands_pilot.json',[])};self.normalization=self.read('normalization.json',{})
  self.scores=pd.read_parquet(self.data/'scores.parquet').set_index('pose_id') if (self.data/'scores.parquet').exists() else pd.DataFrame()
  self.receptor=self.read('receptor_map.json',{});self.structures={};self.ledger={}

  # Keep all 25,900 logical rows visible while payload verification advances.
  for line in self.source_manifest.open():
   r=json.loads(line);r['status']='PAYLOAD_VERIFICATION_PENDING';self.ledger[r['pose_id']]=r
  if (self.data/'status_ledger.jsonl').exists():
   for line in (self.data/'status_ledger.jsonl').open():
    try:r=json.loads(line)
    except json.JSONDecodeError:continue
    self.ledger[r['pose_id']]=r
  self.metric_labels=dict(METRIC_LABELS)
  from .enrichment import configure as configure_enrichment
  configure_enrichment(self)
  from .historical import configure as configure_historical
  configure_historical(self)
  configure(self,self.data)
  self.projection_version=hashlib.sha256((self.data/'receptor_map.json').read_bytes()+(self.data/'enrichment_receptor.json').read_bytes()+VERSION.encode()).hexdigest() if (self.data/'enrichment_receptor.json').exists() else hashlib.sha256((self.data/'receptor_map.json').read_bytes()+VERSION.encode()).hexdigest()
  if (self.data/'enrichment_historical_receptor.json').exists():self.projection_version=hashlib.sha256((self.projection_version+hashlib.sha256((self.data/'enrichment_historical_receptor.json').read_bytes()).hexdigest()).encode()).hexdigest()
  xy=np.array(list(self.receptor.values()),dtype=float);xy=xy[np.isfinite(xy).all(axis=1)]
  self.receptor_projection=dict(basis_sha256=hashlib.sha256((self.data/'receptor_basis.npz').read_bytes()).hexdigest() if (self.data/'receptor_basis.npz').exists() else None,basis_policy='Frozen original fitting population and mean/basis; new receptors transform only; no refit',scope='Finite dashboard receptors, all ligands, raw and relaxed; frozen original shared PCA frame, not refitted; bounds from dashboard receptors only',count=len(xy),bounds=dict(x=[float(xy[:,0].min()),float(xy[:,0].max())],y=[float(xy[:,1].min()),float(xy[:,1].max())]) if len(xy) else None)
  # Expected immutable source bindings are computed once, never inferred from energy rows.
  self.energy_bindings={}
  for pid,source in self.ledger.items():
   if source.get('relaxed_available'):
    try:self.energy_bindings[pid]={'receipt.json':source['receipt_sha256'],**{str(Path(source[field]).relative_to(source['packet_dir'])):source['source_hashes'][source[field]] for field in ('coordinates_path','relaxed_structure_path','system_xml_path')}}
    except (KeyError,ValueError):pass
  from .energy_snapshot import EnergySnapshot
  self.energy_cache=EnergySnapshot(self.data,self.validate_energies)
 def validate_energies(self,f):
  values={}
  for r in f.to_dict('records'):
   r=clean(r);pid=r['pose_id']
   if pid not in self.ledger:continue
   if r.get('energy_status') in ('AVAILABLE','STRAIN_AVAILABLE'):
    hashes=json.loads(r['source_hashes']) if isinstance(r.get('source_hashes'),str) else (r.get('source_hashes') or {})
    expected=self.energy_bindings.get(pid)
    if not expected or any(hashes.get(k)!=v for k,v in expected.items()):
     r.update(energy_status='SOURCE_BINDING_UNAVAILABLE',strain_status='SOURCE_BINDING_UNAVAILABLE',interaction_status='SOURCE_BINDING_UNAVAILABLE',energy_reason='Energy row does not close against this verified coordinate snapshot')
     for col in ('strain_kcal_mol','interaction_kcal_mol','protein_interaction_kcal_mol','heme_interaction_kcal_mol'):r[col]=None
   values[pid]=r
  from .historical import energies as historical_energies
  return historical_energies(self,values)
 def energies(self):return dict(self.energy_cache.read()[0],**self.water_energy)
 def features(self,ligand,state):
  if ligand not in self.ligands:raise KeyError('Unknown ligand')
  if state not in ('raw','relaxed','paired'):raise ValueError('state must be raw, relaxed or paired')
  states=['raw','relaxed'] if state=='paired' else [state];rows=[];coords=[];fe=[]
  pairs=None
  if state=='paired':pairs={r['id'] for r in self.features(ligand,'raw')[0]} & {r['id'] for r in self.features(ligand,'relaxed')[0]}
  for st in states:
   from .historical import feature_folder
   folder=feature_folder(self.data,ligand,st)
   suffix='' if self.full_features else '.pilot'
   r=json.loads((folder/(st+suffix+'.json')).read_text())
   with np.load(folder/(st+suffix+'.npz')) as z:
    coordinates=z['coordinates'];fe_distances=z['fe_distances']
   for i,row in enumerate(r):
    if not included(row['method']):continue
    if pairs is not None and row['id'] not in pairs:continue
    rows.append(row);coords.append(coordinates[i]);fe.append(fe_distances[i])
    if row.get('structure_path'):self.structures[row['structure_key']]=(row['structure_path'],row['structure_sha256'])
  return rows,np.array(coords),np.array(fe)
 def analysis(self,ligand,state,metric,atoms,threshold,subset):
  with COMPUTE:
   if not np.isfinite(threshold) or threshold<0:raise ValueError('threshold must be finite and nonnegative')
   rows,x,fe=self.features(ligand,state);ids=[(r['id']+'::'+r['source_state'] if state=='paired' else r['id']) for r in rows]
   if subset:
    unknown=set(subset)-set(ids)
    if unknown:raise ValueError('pose_ids outside current ligand/state')
    ix=[i for i,p in enumerate(ids) if p in set(subset)];rows=[rows[i] for i in ix];x=x[ix];fe=fe[ix];ids=[ids[i] for i in ix]
   omitted=[]
   if metric in ('contacts','heme'):
    ix=[i for i,r in enumerate(rows) if (r.get('contacts_available',True) if metric=='contacts' else np.isfinite(fe[i]).all())]
    omitted=[pid for i,pid in enumerate(ids) if i not in ix];rows=[rows[i] for i in ix];x=x[ix];fe=fe[ix];ids=[ids[i] for i in ix]
   code_hash=hashlib.sha256((ROOT/'backend/science.py').read_bytes()).hexdigest()
   import scipy,rdkit
   numerical_version=(np.__version__,scipy.__version__,rdkit.__version__,hashlib.sha256((ROOT/'backend/main.py').read_bytes()).hexdigest())
   coordinate_hash=hashlib.sha256(x.tobytes()).hexdigest()
   key=hashlib.sha256(json.dumps([self.snapshot.get('id'),VERSION,code_hash,numerical_version,PROJECTION_VERSION,CLUSTER_ALGORITHM,coordinate_hash,hashlib.sha256(fe.tobytes()).hexdigest(),self.ligands[ligand]['smiles'],ligand,state,metric,atoms,subset,omitted,[(r['id'],r['source_state'],r['structure_sha256'],r.get('contacts'),r.get('contacts_available')) for r in rows]]).encode()).hexdigest()[:24]
   folder=getattr(self,'analysis_folder',self.data/('analysis_cache_'+VERSION));folder.mkdir(exist_ok=True)
   matrix=folder/(key+'.npy');phase=folder/(key+'.npz');precomputed=None
   reused=matrix.exists() and phase.exists()
   if reused:
    try:
     d=np.load(matrix,allow_pickle=False)
     with np.load(phase,allow_pickle=False) as z:
      precomputed={k:z[k] for k in ('xy','linkage')};precomputed['projection']=json.loads(z['projection'].item())
    except (OSError,ValueError,KeyError):reused=False
   if reused:pass
   elif metric=='contacts':d=contact_matrix([r['contacts'] for r in rows])
   elif metric in ('rmsd','selected_atoms'):
    if metric=='selected_atoms' and not atoms:raise ValueError('Select at least one source atom')
    d=rmsd_matrix(x,self.ligands[ligand]['smiles'],atoms if metric=='selected_atoms' else None)
   elif metric=='heme':
    from scipy.spatial.distance import pdist,squareform
    d=rmsd_matrix(fe[:,:,None],self.ligands[ligand]['smiles']) if len(fe) else np.zeros((0,0))
   else:raise ValueError('Unknown metric')
   if not reused:
    from .science import _canonical_pcoa
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform
    n=len(ids)
    xy=_canonical_pcoa(d) if n>2 else np.array([[-d[0,1]/2,0],[d[0,1]/2,0]]) if n==2 else np.zeros((n,2))
    precomputed=dict(xy=xy,linkage=linkage(squareform(d,checks=False),method='average') if n>1 else np.empty((0,4)))
    tmp=folder/(key+f'.{os.getpid()}.npy');np.save(tmp,d);os.replace(tmp,matrix)
   result=analyze_distances(d,ids,threshold,metric,precomputed=precomputed);result['unavailable_pose_ids']=omitted
   if not reused:
    tmp=folder/(key+f'.{os.getpid()}.npz');np.savez(tmp,xy=precomputed['xy'],linkage=precomputed['linkage'],projection=json.dumps(result['projection']));os.replace(tmp,phase)
   os.utime(matrix,None);os.utime(phase,None)
   result['distances_url']='/api/distances/'+key;result['coordinate_state']=state;result['snapshot']=self.snapshot.get('id');result['atom_indices']=list(atoms);result['coordinate_sha256']=coordinate_hash;result['kernel_sha256']=code_hash
   result['analysis_spec']=dict(snapshot=self.snapshot.get('id'),source_snapshot=getattr(self,'source_snapshot',self.snapshot.get('id')),ligand=ligand,state=state,metric=metric,atom_indices=list(atoms),pose_ids=ids,threshold=threshold,kernel_sha256=code_hash,projection_version=PROJECTION_VERSION,cluster_algorithm=CLUSTER_ALGORITHM,population_sha256=hashlib.sha256(json.dumps(ids,separators=(',',':')).encode()).hexdigest())
   result['diagnostics']=dict(distance_kernel_calls=0 if reused else 1,projection_kernel_calls=0 if reused else 1,linkage_kernel_calls=0 if reused else 1,phase_cache_hit=reused,phase_cache_key=key)
   return result
 def meta(self):
  methods=sorted({r['method'] for r in self.ledger.values()});energy=self.energies()
  return dict(receptor_projection=self.receptor_projection,snapshot=self.snapshot.get('id'),default_ligand=next(iter(self.ligands),None),default_state='relaxed',ligands=list(self.ligands.values()),methods=[dict(id=m,label=m,color=PALETTE[i%len(PALETTE)]) for i,m in enumerate(methods)],metrics=[dict(id=k,label=v,unit='native',kind='score',direction='higher',coverage=sum(g['challenge'] for g in self.normalization.get('groups',[]) if g['metric']==k)) for k,v in METRIC_LABELS.items()]+self.extra_metrics+[dict(id=k,label=v,unit='kcal/mol',kind='energy',direction='lower') for k,v in [('strain','Relaxed ligand strain'),('interaction','Protein + heme interaction'),('protein','Protein interaction'),('heme','Heme interaction')]],coverage=dict(self.snapshot,features=self.feature_coverage,feature_scope='frozen dashboard feature population',controls=self.normalization.get('control_count',0),energy_rows=len(energy)),capabilities=dict(contacts=True,typed_hydrogen_bonds=False,pi_stacking=False,non_Fe_heme_fingerprint=False,contacts_protocol='Chemistry-aware atom proximity; polar candidates, not validated H bonds; not ProLIF',selected_atoms=True,shared_receptor_map=True,energy=bool(energy),native_confidence='Inherited parent; exact native semantics',normalization=self.normalization.get('id')))
 def ligand(self,ligand,state):
  info=self.ligands.get(ligand)
  if info is None:raise KeyError('Unknown ligand')
  rows,_,_=self.features(ligand,state);a=self.analysis(ligand,state,'rmsd',(),2.,());energy,energy_publication=self.energy_cache.read();out=[]
  location={pid:i for i,pid in enumerate(a['pose_ids'])}
  for i,r in enumerate(rows):
   pose={k:v for k,v in r.items() if k not in ('contact_details','structure_path')};pid=r['id'];s=self.scores.loc[pid] if pid in self.scores.index else {};e=energy.get(pid,{}) if r['source_state']=='relaxed' else {}
   if r['source_state']=='relaxed' and pid in self.water_energy:e=self.water_energy[pid]
   pose.update(scores={k:s.get(k) for k in self.metric_labels},percentiles={k:s.get('pct__'+k) for k in self.metric_labels},energies={k:e.get(col) for k,col in [('strain','strain_kcal_mol'),('interaction','interaction_kcal_mol'),('protein','protein_interaction_kcal_mol'),('heme','heme_interaction_kcal_mol')]},energy_status=e.get('energy_status','unavailable_raw_state' if r['source_state']=='raw' else 'not_yet_scored'),energy_reason=e.get('energy_reason'),pose_xy=None,receptor_xy=self.receptor.get(r['source_state']+'|'+pid),cluster=None,native_confidence_label='Native parent prediction confidence (not recomputed)' if r['source_state']=='relaxed' else 'Native prediction confidence')
   from .enrichment import decorate
   decorate(self,pose,pid)
   pose['strain_status']=e.get('strain_status',e.get('energy_status','not_yet_scored'));pose['interaction_status']=e.get('interaction_status',e.get('energy_status','not_yet_scored'))
   pose['energy_reference_version']=e.get('reference_version');pose['energy_protocol_id']=e.get('protocol_id')
   pose['energy_reference_key']=e.get('reference_key');pose['protein_state_key']=e.get('protein_state_key');pose['energy_reference_search_warning']=e.get('reference_search_warning')
   view_id=pid+'::'+r['source_state'] if state=='paired' else pid
   ai=location.get(view_id);pose['pose_xy']=a['xy'][ai] if ai is not None else None;pose['cluster']=a['labels'][ai] if ai is not None else None
   if state=='paired':pose['parent_id']=pid;pose['id']=view_id
   out.append(pose)
  displayed={p['id'].split('::')[0] for p in out}
  unavailable=[dict(id=r['pose_id'],status=r['status'],reason=r.get('verification_error') or r.get('reason')) for r in self.ledger.values() if r['ligand_id']==ligand and r['pose_id'] not in displayed]
  return clean(dict(analysis_spec=a['analysis_spec'],energy_publication=energy_publication,ligand=ligand,arm=info.get('arm','dry'),ligand_label=info.get('label',ligand),parent_ligand_id=info.get('parent_ligand_id'),state=state,snapshot=self.snapshot.get('id'),smiles=info['smiles'],svg_url=info['svg_url'],atoms=info['atoms'],moieties=info['moieties'],poses=out,coverage=dict(total=info['n_total'],shown=len(out),unavailable=len(unavailable),unavailable_poses=unavailable),clusters=a['clusters'],projection=a['projection'],normalization=dict(id=self.normalization.get('id'),scope=self.normalization.get('scope'))))
 def versions(self):
  return dict(snapshot=self.snapshot.get('id'),projection=self.projection_version,normalization=self.normalization.get('id'),energy_publication=self.energy_cache.read()[1])
 def background(self,ligand,state,hover_metric=None):
  rows,_,_=self.features(ligand,state);out=[]
  for r in rows:
   pid=r['id'];xy=self.receptor.get(r['source_state']+'|'+pid)
   if xy is None:continue
   s=self.scores.loc[pid] if pid in self.scores.index else {}
   pose=dict(id=pid+'::'+r['source_state'] if state=='paired' else pid,method=r['method'],seed=r.get('seed'),sample=r.get('sample'),model=r.get('model'),source_state=r['source_state'],receptor_xy=xy,scores={},percentiles={})
   from .enrichment import decorate
   decorate(self,pose,pid)
   if hover_metric and hover_metric not in ('method','cluster'):
    pose['scores']={hover_metric:pose['scores'].get(hover_metric,s.get(hover_metric))}
    pose['percentiles']={hover_metric:pose['percentiles'].get(hover_metric,s.get('pct__'+hover_metric))}
   else:pose.pop('scores');pose.pop('percentiles')
   # These are the native-confidence fields read by the background hover.
   for k in list(pose):
    if k.startswith('native_confidence_') and k!='native_confidence_percentile':pose.pop(k)
   pose.pop('arm',None);pose.pop('source_pose_id',None)
   if hover_metric in ('ser119_contact','arg105_contact'):
    pose.update(structure_sha256=r.get('structure_sha256'),available=r.get('available',True))
    if getattr(self,'_hbond_ligand',None)!=ligand:
     folder=ROOT/'static/residue-hbonds';manifest_path=folder/'manifest.json'
     if manifest_path.exists():
      manifest=json.loads(manifest_path.read_text());shard=manifest['shards'][ligand];path=folder/shard['file']
      if path.parent!=folder:raise ValueError('Invalid Hbond shard path')
      raw=path.read_bytes()
      if hashlib.sha256(raw).hexdigest()!=shard['sha256']:raise ValueError('Hbond shard hash mismatch')
      payload=json.loads(raw)
      if payload['ligand']!=ligand or any(payload.get(k)!=manifest.get(k) for k in ('source_snapshot','source_scientific_revision','source_index_sha256')):raise ValueError('Hbond shard source mismatch')
      from .historical import feature_folder
      for source_state,identity in shard['feature_sources'].items():
       source=feature_folder(self.data,ligand,source_state)/(source_state+'.json')
       if hashlib.sha256(source.read_bytes()).hexdigest()!=identity['sha256']:raise ValueError('Hbond feature source mismatch')
      self._hbond_rows=payload['poses']
     else:self._hbond_rows=json.loads((folder/'index.json').read_text())['poses']
     self._hbond_ligand=ligand
    row=self._hbond_rows.get(r['source_state']+'|'+pid,{})
    value=row.get(hover_metric) if pose['available'] and pose['structure_sha256'] and row.get('structure_sha256')==pose['structure_sha256'] else None
    pose['contact_metrics']={hover_metric:value if value in (0,1) else None}
   out.append(pose)
  return clean(dict(ligand=ligand,state=state,snapshot=self.snapshot.get('id'),projection_version=self.projection_version,poses=out))
 def pose(self,pid,state):
  pid=pid.split('::')[0];row=self.ledger.get(pid)
  if row is None:raise KeyError('Unknown pose')
  for r in self.features(row['ligand_id'],state)[0]:
   if r['id']==pid:return clean(dict(r,provenance=row,native_scores=self.scores.loc[pid].to_dict() if pid in self.scores.index else None,energy=self.energies().get(pid) if state=='relaxed' else None))
  raise KeyError('Coordinates unavailable in this state')
 def structure(self,key):
  if key not in self.structures:raise KeyError('Unknown structure key')
  p,h=self.structures[key];p=Path(p)
  from analysis.features_freeze import sha,safe
  assert sha(p)==h,'Source structure changed'
  return p.read_bytes()
STORE=None
def get_store():
 global STORE
 if STORE is None:STORE=Store()
 return STORE

def route(path,query=None,body=None):
 s=get_store();query=query or {};state=query.get('state','relaxed')
 if path=='/api/meta':return s.meta()
 if path=='/api/versions':return s.versions()
 if path=='/api/receptor-background':return s.background(query.get('ligand'),state,query.get('hover_metric'))
 if path=='/api/health':return dict(status='ok',snapshot=s.snapshot,ligands=len(s.ligands),feature_build=s.feature_coverage)
 if path=='/api/ligand':return s.ligand(query.get('ligand',next(iter(s.ligands))),state)
 if path=='/api/pose':return s.pose(query.get('id',''),state)
 if path=='/api/analyze':
  b=body or {};atoms=b.get('atoms',[])
  if not isinstance(atoms,list) or any(type(a) is not int for a in atoms):raise ValueError('atoms must be integer source indices')
  return s.analysis(b.get('ligand'),b.get('state','relaxed'),b.get('metric','contacts'),tuple(sorted(set(atoms))),float(b.get('threshold',.55)),tuple(b.get('pose_ids') or []))
 if path.startswith('/api/structure/'):return s.structure(path.rsplit('/',1)[1]),'chemical/x-pdb'
 if path.startswith('/api/molecule/'):
  name=path.rsplit('/',1)[1];ligand=name.removesuffix('.svg')
  if ligand not in s.ligands or name!=ligand+'.svg':raise KeyError('Unknown molecule')
  return (DATA/'molecules'/name).read_bytes(),'image/svg+xml'
 if path.startswith('/api/distances/'):
  key=path.rsplit('/',1)[1]
  if len(key)!=24 or any(c not in '0123456789abcdef' for c in key):raise KeyError('Unknown matrix')
  p=DATA/('analysis_cache_'+VERSION)/(key+'.npy')
  if not p.exists():raise KeyError('Unknown matrix')
  return p.read_bytes(),'application/octet-stream'
 if path=='/api/schema':return read('api_schema.json',{})
 if path=='/':p=ROOT/'static/index.html'
 elif path.startswith('/static/'):p=ROOT/path.lstrip('/')
 else:raise KeyError('Unknown route')
 if '..' in Path(path).parts or not p.resolve().is_relative_to((ROOT/'static').resolve()) or not p.is_file():raise KeyError('Unknown asset')
 return p.read_bytes(),mimetypes.guess_type(p.name)[0] or 'application/octet-stream'

if __name__=='__main__':
 import uvicorn
 parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=8828);parser.add_argument('--host',default='0.0.0.0');args=parser.parse_args()
 uvicorn.run('backend.asgi:app',host=args.host,port=args.port,workers=1,limit_concurrency=32,backlog=32,timeout_keep_alive=3,timeout_graceful_shutdown=8,h11_max_incomplete_event_size=16384)
