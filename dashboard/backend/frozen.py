"""Original Store methods with exported immutable, per-ligand inputs."""
import json,gzip,types
from pathlib import Path
import pandas as pd
from .main import Store,ROOT,DATA
class FrozenStore(Store):
 def __init__(self,ligand=None):
  self.data=DATA
  self.read=lambda name,default=None:json.loads((DATA/name).read_text()) if (DATA/name).exists() else default
  props=json.loads((DATA/'frozen/properties.json').read_text())
  for k,v in props.items():setattr(self,k,v)
  self.source_snapshot=json.loads((ROOT/'manifest.json').read_text())['snapshot']
  if ligand is None:ligand=next(iter(self.ligands))
  if ligand not in self.ligands:raise KeyError('Unknown ligand')
  self.full_features=True;self.structures={};self.water_energy={};self.historical_energy={}
  self.scores=pd.read_parquet(DATA/'frozen'/(ligand+'.scores.parquet'))
  energy=json.loads((DATA/'frozen'/(ligand+'.energies.json')).read_text())
  self.energy_cache=types.SimpleNamespace(read=lambda:(energy,self.energy_status))
  self.analysis_folder=ROOT/'cache/analysis';self.analysis_folder.mkdir(parents=True,exist_ok=True)
 def analysis(self,*args,**kwargs):
  from .revealed import decorate_analysis
  return decorate_analysis(self,super().analysis(*args,**kwargs))
 def ligand(self,*args,**kwargs):
  result=super().ligand(*args,**kwargs)
  if hasattr(self,'_revealed_reference'):result['revealed_reference']=self._revealed_reference
  return result
 def pose(self,pid,state):
  result=super().pose(pid,state);lid=self.ledger[pid.split('::')[0]]['ligand_id']
  with gzip.open(DATA/'frozen'/(lid+'.provenance.json.gz'),'rt') as f:result['provenance']=json.load(f)[pid.split('::')[0]]
  with gzip.open(DATA/'features'/lid/(result['source_state']+'.details.json.gz'),'rt') as f:result['contact_details']=json.load(f).get(pid.split('::')[0],[])
  return result
