"""Disposable original-kernel worker with validated immutable inputs."""
import sys,json,math,resource,time
from .frozen import FrozenStore
from .main import clean

def run(job):
 op=job['operation'];q=job.get('query',{});b=job.get('body') or {};lid=b.get('ligand') if op=='analyze' else q.get('ligand')
 if op=='pose':
  from .main import DATA
  props=json.loads((DATA/'frozen/properties.json').read_text());row=props['ledger'].get(q.get('id','').split('::')[0]);lid=row['ligand_id'] if row else None
 s=FrozenStore(lid);state=q.get('state','relaxed')
 import os
 receptor31=os.environ.get('CYP_RECEPTOR31_DIR')
 if receptor31:
  from .receptor31_adapter import configure
  configure(s,receptor31)
 from .truth_metrics import configure as configure_truth_metrics
 configure_truth_metrics(s)
 if op=='ligand':
  result=s.ligand(lid or next(iter(s.ligands)),state)
  fmt=q.get('format')
  if fmt not in (None,'legacy','columns-v1'):raise ValueError('Unknown ligand format')
  if fmt=='columns-v1':
   from .compact import encode
   result=encode(result)
 elif op=='receptor-background':result=s.background(lid,state,q.get('hover_metric'))
 elif op=='pose':result=s.pose(q.get('id',''),state)
 elif op=='analyze':
  atoms=b.get('atoms',[]);subset=b.get('pose_ids') or [];threshold=float(b.get('threshold',2.))
  if not isinstance(atoms,list) or len(atoms)>512 or any(type(a) is not int or a<0 for a in atoms):raise ValueError('Invalid atom indices')
  if not isinstance(subset,list) or len(subset)>5000 or any(not isinstance(p,str) or len(p)>256 for p in subset):raise ValueError('Invalid pose IDs')
  if not math.isfinite(threshold) or threshold<=0:raise ValueError('Threshold must be finite and positive')
  result=s.analysis(lid,b.get('state','relaxed'),b.get('metric','rmsd'),tuple(sorted(set(atoms))),threshold,tuple(subset))
  # Preserve original top-level ordered IDs without repeating them in metadata.
  spec=dict(result['analysis_spec']);spec.pop('pose_ids');spec['pose_ids_source']='pose_ids';result={**result,'analysis_spec':spec}
 else:raise KeyError('Unknown operation')
 return clean(result)
if __name__=='__main__':
 start=time.monotonic()
 try:result={'status':200,'result':run(json.load(sys.stdin))}
 except (ValueError,TypeError,IndexError) as e:result={'status':400,'detail':str(e)}
 except (KeyError,FileNotFoundError) as e:result={'status':404,'detail':'Unknown resource'}
 result['worker_peak_rss_kib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;result['seconds']=time.monotonic()-start
 sys.stdout.write(json.dumps(result,separators=(',',':'),allow_nan=False))
