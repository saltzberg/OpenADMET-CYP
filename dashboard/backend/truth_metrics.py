"""Sealed official geometry metrics, strictly bound to this frozen display universe.

BiSymRMSD (UI) maps explicitly to official source key BiSyRMSD. No scoring,
penalty, percentile conversion, or raw/relaxed fallback occurs here.
"""
import collections
import hashlib
import json
from pathlib import Path
import types

# Filled from the deterministic packet build; manifest itself is also sealed.
MANIFEST_SHA256='cf24bbb7fb39e53e9c6ca991ae43d2d58569f6edd5df524280c8e12252d4c84f'
FIELDS=('LDDT-PLI','BiSyRMSD')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def verify(root):
 root=Path(root);folder=root/'data/truth-metrics';manifest=folder/'manifest.json'
 if sha(manifest)!=MANIFEST_SHA256:raise ValueError('Truth metric manifest seal mismatch')
 m=json.loads(manifest.read_text())
 if json.loads((root/'manifest.json').read_text())['snapshot']!=m['source_snapshot']:raise ValueError('Truth metric scientific snapshot mismatch')
 for name,digest in m['source_guards'].items():
  if sha(root/name)!=digest:raise ValueError('Truth metric frozen source guard mismatch: '+name)
 for name,digest in m['files'].items():
  if sha(folder/name)!=digest:raise ValueError('Truth metric packet hash mismatch: '+name)
 if {p.name for p in folder.iterdir() if p.is_file()}!={'manifest.json',*m['files']}:raise ValueError('Unexpected truth metric packet file')
 return m

def revision_paths(root):
 """Include every derivative file in runtime SID and active revision guards."""
 m=verify(root);folder=Path(root)/'data/truth-metrics'
 return [folder/'manifest.json',*[folder/n for n in sorted(m['files'])]]

class Packet:
 def __init__(self,root):
  self.root=Path(root);self.manifest=verify(root);self.loaded={}
 def rows(self,ligand):
  if ligand not in self.loaded:
   name=ligand+'.json'
   if name not in self.manifest['files']:raise ValueError('Truth metric unknown ligand')
   items=json.loads((self.root/'data/truth-metrics'/name).read_text());rows={}
   for r in items:
    if r['reference_revision']!=self.manifest['reference_revision'] or r['reference_sha256']!=self.manifest['references'][ligand]:raise ValueError('Truth reference identity mismatch')
    key=(r['source_state'],r['pose_id'])
    if key in rows or r['ligand_id']!=ligand:raise ValueError('Truth metric duplicate or wrong ligand')
    rows[key]=r
   self.loaded[ligand]=rows
  return self.loaded[ligand]
 def decorate(self,pose,ligand):
  state=pose.get('source_state');view_id=pose.get('id','');pid=view_id.removesuffix('::'+str(state));r=self.rows(ligand).get((state,pid))
  valid=r is not None and pose.get('structure_sha256')==r['structure_sha256'] and pose.get('method')==r['method'] and r['ligand_id']==ligand and r['reference_revision']==self.manifest['reference_revision'] and pose.get('available',True) is not False and r['reference_sha256']==self.manifest['references'].get(ligand)
  if not valid:
   pose['truth_metrics']={**dict.fromkeys(FIELDS),'status':'unavailable','missing_reason':'display_coordinate_binding_mismatch'}
  else:
   pose['truth_metrics']={k:v for k,v in r.items() if k not in {'smiles','structure_url','coordinate_path','model_path','displayed_model_path','ref_path'}}
   pose['truth_metrics']['status']='scored' if all(r[k] is not None for k in FIELDS) else 'unavailable'
  return pose

def configure(store,root=None):
 if getattr(store,'_truth_metrics_configured',False):raise ValueError('Truth metrics already configured')
 packet=Packet(root or Path(__file__).resolve().parents[1]);ligand_original=store.ligand;pose_original=store.pose
 def ligand(self,lid,state):
  result=ligand_original(lid,state)
  for pose in result['poses']:packet.decorate(pose,lid)
  result['truth_metric_coverage']=dict(collections.Counter(p['truth_metrics']['missing_reason'] or 'scored' for p in result['poses']))
  return result
 def pose(self,pid,state):
  result=pose_original(pid,state);lid=self.ledger[pid.split('::')[0]]['ligand_id'];return packet.decorate(result,lid)
 store.ligand=types.MethodType(ligand,store);store.pose=types.MethodType(pose,store);store._truth_metrics_configured=True
