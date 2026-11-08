"""Prefix-scoped read-only API. No campaign reads, mutable routes or arbitrary files."""
import asyncio,json,hashlib,re,mimetypes,os,sys,time
from pathlib import Path
from fastapi import FastAPI,Request
from fastapi.responses import Response
from starlette.exceptions import HTTPException
from .workers import Workers,Busy,Deadline,WorkerFailed
from .runtime_cache import DiskCache,Shared,canonical_key
ROOT=Path(__file__).resolve().parents[1];PREFIX='/articles/admet-cyp-structure/cofolding/'
manifest=json.loads((ROOT/'manifest.json').read_text());SOURCE_SNAPSHOT=manifest['snapshot']
structures=json.loads((ROOT/'structures.json').read_text());assets={v['path'] for v in structures.values()}
static={str(p.relative_to(ROOT/'static')):p for p in (ROOT/'static').rglob('*') if p.is_file()}
viewer_references=json.loads((ROOT/'static/revealed/index.json').read_text())
# Hash a bounded, explicit application tree once at startup. Scientific source
# identity stays separate; neither the scientific manifest nor inputs change.
reference_paths=sorted((ROOT/'data/revealed-reference').glob('*.json'))
revision_paths=sorted(set([ROOT/'manifest.json',ROOT/'structures.json',ROOT/'data/api_schema.json']+reference_paths+list((ROOT/'backend').glob('*.py'))+list(static.values())))
RECEPTOR31_DIR=Path(os.environ['CYP_RECEPTOR31_DIR']) if os.environ.get('CYP_RECEPTOR31_DIR') else None
RECEPTOR31_INFO=None
if RECEPTOR31_DIR is not None:
 if not RECEPTOR31_DIR.resolve().is_relative_to(ROOT.resolve()):raise ValueError('Stage receptor31 data inside the application root')
 from .receptor31_adapter import scientific_revision
 receptor31_revision=scientific_revision(RECEPTOR31_DIR)
 RECEPTOR31_INFO=dict(json.loads((RECEPTOR31_DIR/'manifest.json').read_text()),revision=receptor31_revision,scientific_source_snapshot=SOURCE_SNAPSHOT)
 revision_paths += [RECEPTOR31_DIR/n for n in ('manifest.json','receptor_basis.npz','receptor_map.json','revealed_map.json','artifact-hashes.json')]
from .truth_metrics import revision_paths as truth_revision_paths
revision_paths += truth_revision_paths(ROOT)
REVISION={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in revision_paths if p.is_file()}
REVISION_STATS={name:((ROOT/name).stat().st_mtime_ns,(ROOT/name).stat().st_size,(ROOT/name).stat().st_ino) for name in REVISION}
SID=hashlib.sha256(json.dumps([SOURCE_SNAPSHOT,REVISION],sort_keys=True,separators=(',',':')).encode()).hexdigest();BASE=PREFIX+'api/v1/'+SID
DIAGNOSTICS=os.environ.get('CYP_DIAGNOSTICS')=='1'
workers=Workers(concurrency=1,queued=2,deadline=45,diagnostics=DIAGNOSTICS)
CACHE=ROOT/'cache';CACHE.mkdir(exist_ok=True);MAX_CACHE=128*1024**2
cache=DiskCache(CACHE,MAX_CACHE,diagnostics=DIAGNOSTICS);shared=Shared(cache)
COORDS=Path(os.environ.get('CYP_COORDINATE_ROOT',str(ROOT/'smoke-coordinates')))
app=FastAPI(docs_url=None,redoc_url=None,openapi_url=None)
app.router.add_event_handler('shutdown',workers.close)

def encode(x):return json.dumps(x,separators=(',',':'),allow_nan=False).encode()
def sanitize(x):
 if isinstance(x,dict):return {(('private-source-'+hashlib.sha256(k.encode()).hexdigest()[:16]) if isinstance(k,str) and k.startswith(('/home/','/archive/','/mnt/','/data')) else k):sanitize(v) for k,v in x.items() if k not in {'structure_path','packet_dir','confidence_path','traceback'}}
 if isinstance(x,list):return [sanitize(v) for v in x]
 if isinstance(x,str):
  if x.startswith('/api/'):return BASE+x[4:]
  if x.startswith(('/home/','/archive/','/mnt/','/data','/opt/')):return '[private source path]'
  # Tracebacks/errors can embed a path rather than start with one.
  if re.search(r'/(?:home|archive|mnt|data[0-9]*)/',x):return '[private source provenance]'
 return x

def reply(value,status=200,ctype='application/json',pinned=False,request=None,revalidate=False,headers=None):
 payload=value if isinstance(value,bytes) else encode(value)
 policy=('public, max-age=31536000, immutable' if pinned else 'no-cache') if status==200 and request is not None and request.method=='GET' and (pinned or revalidate) else 'no-store'
 h={'Cache-Control':policy,'X-Content-Type-Options':'nosniff','X-Snapshot':SID,**(headers or {})}
 if policy!='no-store':
  assert request is not None
  etag='"'+hashlib.sha256(payload).hexdigest()+'"';h['ETag']=etag
  # If-None-Match uses weak comparison for GET, including a list or wildcard.
  candidates=[x.strip().removeprefix('W/') for x in request.headers.get('if-none-match','').split(',')]
  if '*' in candidates or etag in candidates:return Response(status_code=304,headers=h)
 return Response(payload,status_code=status,media_type=ctype,headers=h)

@app.exception_handler(HTTPException)
async def http_error(request,exception):
 return reply({'detail':exception.detail},exception.status_code,headers=exception.headers)

def trim():
 return cache.trim()

def active_revision():
 # Detect in-place development/deployment writes without rehashing vendor/data
 # bundles on each request. Published workers must use the startup revision.
 try:return {str(p.relative_to(ROOT)) for p in (ROOT/'data/truth-metrics').iterdir() if p.is_file()}=={n for n in REVISION if n.startswith('data/truth-metrics/')} and all(((ROOT/name).stat().st_mtime_ns,(ROOT/name).stat().st_size,(ROOT/name).stat().st_ino)==stamp for name,stamp in REVISION_STATS.items()) and {str(p.relative_to(ROOT)) for p in (ROOT/'backend').glob('*.py')}=={n for n in REVISION if n.startswith('backend/')} and {str(p.relative_to(ROOT)) for p in (ROOT/'data/revealed-reference').glob('*.json')}=={n for n in REVISION if n.startswith('data/revealed-reference/')}
 except FileNotFoundError:return False

def bind_revision(result):
 result=sanitize(result)
 if not isinstance(result,dict):raise WorkerFailed('Expected worker result object')
 source=result.get('snapshot')
 for field in ('analysis_spec','columns','revealed_reference'):
  nested=result.get(field)
  if isinstance(nested,dict):
   if nested.get('snapshot')!=source:raise WorkerFailed('Mixed scientific response revision')
   nested['source_snapshot']=nested.get('source_snapshot',source)
   nested['snapshot']=SID
   if field=='revealed_reference':nested['viewer_pose']=viewer_references.get(nested.get('ligand'))
 result.update(snapshot=SID,source_snapshot=SOURCE_SNAPSHOT,energy_revision=manifest['source_snapshot']['energy_revision'])
 return result

async def compute(op,q,body,key):
 timing={};job=encode({'operation':op,'query':q,'body':body})
 if not active_revision():raise WorkerFailed('Application revision changed; restart required')
 with cache.stage(seed=False) as (folder,existing):
  raw=await workers.execute(job,command=[sys.executable,'-m','backend.runtime_cache'],background=op=='receptor-background',timing=timing,env={'CYP_RUNTIME_ANALYSIS_STAGE':str(folder),'CYP_RUNTIME_ANALYSIS_CACHE':str(cache.analysis)})
  if (folder/'.seeds').is_file():existing={name:tuple(value) for name,value in json.loads((folder/'.seeds').read_text()).items()}
  start=time.monotonic()
  try:result=json.loads(raw)
  except (ValueError,TypeError):raise WorkerFailed('Invalid worker response')
  timing['decode']=time.monotonic()-start
  if not active_revision():raise WorkerFailed('Application revision changed; restart required')
  if result['status']!=200:return None,result['status'],sanitize({'detail':result['detail']}),timing
  start=time.monotonic();payload=encode(bind_revision(result['result']));timing['serialization']=time.monotonic()-start
  cache.count('serialization_seconds',timing['serialization'])
  # Pin the current job's scientific files during publication and incoming
  # response reservation. Stage hard links keep queued readers safe on eviction.
  distance_key=str(result['result'].get('distances_url','')).rsplit('/',1)[-1]
  pins={cache.analysis/p.name for p in folder.iterdir() if p.suffix in {'.json','.npy','.npz'} and (existing.get(p.name)!=(p.stat().st_ino,p.stat().st_size) or p.stem==distance_key)}
  cache.protected.update(pins)
  try:
   cache.promote(folder,existing)
   if len(payload)<=32*1024**2:cache.put(CACHE/(key+'.json'),payload)
  finally:cache.protected.difference_update(pins)
  return payload,200,None,timing
@app.api_route('/{path:path}',methods=['GET','POST','PUT','DELETE','PATCH','OPTIONS','HEAD','TRACE'])
async def endpoint(path:str,request:Request):
 try:
  if len(request.url.path)>2048 or len(request.url.query)>8192:return reply({'detail':'URL too long'},414)
  if request.method not in {'GET','POST'}:return reply({'detail':'Method not allowed'},405)
  full='/'+path
  if not full.startswith(PREFIX):return reply({'detail':'Unknown route'},404)
  p=full[len(PREFIX):]
  if p=='api/bootstrap':
   if request.method!='GET':return reply({'detail':'Method not allowed'},405)
   if not active_revision():raise WorkerFailed('Application revision changed; restart required')
   return reply({'snapshot':SID,'source_snapshot':SOURCE_SNAPSHOT,'scientific_source':manifest['source_snapshot'],'api_base':BASE,'energy_revision':manifest['source_snapshot']['energy_revision'],'analysis_kernel_sha256':REVISION['backend/science.py'],'projection_version':'pcoa-canonical-eigenspaces-v2','cluster_algorithm':'average-linkage-original-distances-v1','asset_policy':'revalidate','asset_version':SID,'old_tab_policy':'Old runtime snapshot URLs return 404; reload required'})
  if p.startswith('api/'):
   prefix='api/v1/'+SID+'/'
   if not p.startswith(prefix):return reply({'detail':'Unknown runtime snapshot; reload required'},404)
   op=p[len(prefix):]
   if (request.method=='POST')!=(op=='analyze'):return reply({'detail':'Method not allowed'},405)
   pairs=list(request.query_params.multi_items());q=dict(pairs);body=None
   if op=='analyze':
    raw=bytearray()
    async with asyncio.timeout(5):
     async for chunk in request.stream():
      raw.extend(chunk)
      if len(raw)>65536:return reply({'detail':'Request exceeds 64 KiB'},413)
    body=json.loads(raw)
    if not isinstance(body,dict):raise ValueError('Expected JSON object')
   if op=='health':return reply({'status':'ok','snapshot':SID})
   if op=='diagnostics' and DIAGNOSTICS:
    return reply({'workers':workers.counters,'cache':cache.counters,'cache_bytes':cache.size(),'cache_limit':cache.limit,'admitted':workers.admitted,'active':int(workers.active is not None),'foreground_queued':len(workers.foreground),'background_admitted':workers.background_admitted,'inflight':len(shared.entries),'inflight_waiters':sum(e['waiters'] for e in shared.entries.values())})
   if not active_revision():raise WorkerFailed('Application revision changed; restart required')
   if op in {'meta','versions','schema'}:
    target=ROOT/'data'/('api_schema.json' if op=='schema' else 'frozen/'+op+'.json');result=json.loads(target.read_text())
    if RECEPTOR31_INFO is not None:
     if op=='meta':result=dict(result,receptor_projection=RECEPTOR31_INFO)
     elif op=='versions':result=dict(result,projection=RECEPTOR31_INFO['projection_version']+':'+RECEPTOR31_INFO['revision'])
   elif op.startswith('structure/'):
    key=op.split('/',1)[1]
    if key not in structures:raise KeyError()
    return Response(status_code=307,headers={'Location':PREFIX+structures[key]['path'],'Cache-Control':'no-store'})
   elif op.startswith('molecule/'):
    name=op.split('/',1)[1]
    if not re.fullmatch(r'OCNT-[0-9]+\.svg',name):raise KeyError()
    target=ROOT/'data/molecules'/name
    if not target.is_file():raise KeyError()
    return reply(target.read_bytes(),ctype='image/svg+xml',pinned=True,request=request)
   elif op.startswith('distances/'):
    key=op.split('/',1)[1]
    if not re.fullmatch('[0-9a-f]{24}',key):raise KeyError()
    target=CACHE/'analysis'/(key+'.npy')
    payload=cache.read(target)
    if payload is None:raise KeyError()
    return reply(payload,ctype='application/octet-stream',pinned=True,request=request)
   elif op in {'ligand','pose','receptor-background','analyze'}:
    if not active_revision():raise WorkerFailed('Application revision changed; restart required')
    key=canonical_key(SID,op,pairs,body);target=CACHE/(key+'.json');payload=cache.read(target)
    if payload is not None and op=='analyze':
     # An analysis response must not strand a download URL after its matrix
     # was evicted. Repeating the same POST then rebuilds the missing unit.
     distance_key=str(json.loads(payload).get('distances_url','')).rsplit('/',1)[-1]
     if not re.fullmatch('[0-9a-f]{24}',distance_key) or not (cache.analysis/(distance_key+'.npy')).is_file():
      target.unlink(missing_ok=True);payload=None;cache.count('hits',-1);cache.count('misses')
    origin='hit';timing={}
    if payload is None:
     origin='shared' if key in shared.entries else 'miss'
     payload,status,error,timing=await shared.run(key,lambda:compute(op,q,body,key),request=request)
     if status!=200:return reply(error,status)
    headers={'X-Origin-Cache':origin,'Server-Timing':', '.join(f'{k};dur={v*1000:.3f}' for k,v in timing.items())} if DIAGNOSTICS else None
    return reply(payload,pinned=request.method=='GET',request=request,headers=headers)
   else:raise KeyError()
   return reply(bind_revision(result),pinned=True,request=request)
  if request.method!='GET':return reply({'detail':'Method not allowed'},405)
  if p.startswith('cyp-main/'):
   if p not in assets:raise KeyError()
   target=COORDS/p
   if not target.is_file():return reply({'detail':'Coordinate not installed in this candidate subset'},404)
   # gzip bytes as an opaque object: NO Content-Encoding, client decodes once.
   return reply(target.read_bytes(),ctype='application/gzip',revalidate=True,request=request)
  name='index.html' if not p else p.removeprefix('static/')
  if name not in static:raise KeyError()
  if not active_revision():raise WorkerFailed('Application revision changed; restart required')
  return reply(static[name].read_bytes(),ctype=mimetypes.guess_type(name)[0] or 'application/octet-stream',revalidate=True,request=request)
 except (KeyError,FileNotFoundError):return reply({'detail':'Unknown resource'},404)
 except (ValueError,TypeError):return reply({'detail':'Invalid request'},400)
 except Busy:return reply({'detail':'Analysis queue full; retry later'},429,headers={'Retry-After':'1'})
 except (Deadline,TimeoutError):return reply({'detail':'Analysis exceeded 45 second execution/queue deadline'},504)
 except WorkerFailed:return reply({'detail':'Scientific worker failed or application revision changed; retry after reload'},503,headers={'Retry-After':'1'})
 except Exception:return reply({'detail':'Origin request failed; retry later'},503,headers={'Retry-After':'1'})
