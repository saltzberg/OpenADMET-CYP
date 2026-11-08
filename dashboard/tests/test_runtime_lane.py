"""Run from the workspace with the existing service Python; no network required.

All ASGI roots, caches, source variants and fixture commands are isolated copies.
Real scientific requests use the active source worker and frozen input files.
"""
import asyncio
import hashlib
import gzip
import importlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
import types
import unittest
from unittest.mock import patch

import httpx

APP=Path(os.environ.get('CYP_TEST_APP', Path(__file__).resolve().parents[1])).resolve()
sys.path.insert(0,str(APP))
from backend.workers import Workers,Busy,Deadline,WorkerFailed
from backend.runtime_cache import DiskCache,Shared,canonical_key
SCRATCH=Path(os.environ['TMPDIR'])


def command(delay=.05,fail=False,value='ok'):
 return [sys.executable,'-c',f'import sys,time; sys.stdin.buffer.read(); time.sleep({delay!r}); sys.stdout.write({value!r}); sys.exit({int(fail)})']


class WorkerTests(unittest.IsolatedAsyncioTestCase):
 async def test_priority_and_background_admission(self):
  w=Workers(diagnostics=True);order=[]
  async def run(name,bg=False,delay=.05):
   await w.execute(b'',command(delay,value=name),background=bg);order.append(name)
  first=asyncio.create_task(run('active',delay=.2));await asyncio.sleep(.03)
  bg=asyncio.create_task(run('background',True));await asyncio.sleep(.01)
  with self.assertRaises(Busy):await w.execute(b'',command(),background=True)
  a=asyncio.create_task(run('foreground-a'));await asyncio.sleep(.01)
  self.assertEqual(w.admitted,3)
  with self.assertRaises(Busy):await w.execute(b'',command())
  await first
  b=asyncio.create_task(run('foreground-b'));await asyncio.sleep(.01)
  self.assertLessEqual(w.admitted,3)
  await asyncio.gather(first,bg,a,b)
  self.assertEqual(order,['active','foreground-a','foreground-b','background'])
  self.assertEqual(w.admitted,0);self.assertFalse(w.pids)
  self.assertEqual(w.counters['launches'],4)
  print('PRIORITY',json.dumps({'order':order,'counters':w.counters}))
 async def test_active_background_allows_two_foreground(self):
  w=Workers();bg=asyncio.create_task(w.execute(b'',command(.15),background=True));await asyncio.sleep(.02)
  a=asyncio.create_task(w.execute(b'',command()));b=asyncio.create_task(w.execute(b'',command()));await asyncio.sleep(.02)
  with self.assertRaises(Busy):await w.execute(b'',command())
  await asyncio.gather(bg,a,b);self.assertEqual(w.admitted,0)
 async def test_errors_timeout_queued_cancel_and_retry(self):
  w=Workers(deadline=.12,diagnostics=True)
  with self.assertRaises(WorkerFailed):await w.execute(b'',command(.01,True))
  with self.assertRaises(Deadline):await w.execute(b'',command(1))
  self.assertFalse(w.pids);self.assertEqual(w.admitted,0)
  self.assertEqual(await w.execute(b'',command(.01)),b'ok')
  active=asyncio.create_task(w.execute(b'',command(1)));await asyncio.sleep(.02)
  queued=asyncio.create_task(w.execute(b'',command()));await asyncio.sleep(.01);queued.cancel()
  with self.assertRaises(asyncio.CancelledError):await queued
  active.cancel()
  with self.assertRaises(asyncio.CancelledError):await active
  self.assertFalse(w.pids);self.assertEqual(w.admitted,0);self.assertFalse(w.foreground)
 async def test_stdout_limit_remains_32_mib(self):
  w=Workers();cmd=[sys.executable,'-c','import sys; sys.stdout.buffer.write(b"x"*(32*1024**2+1))']
  with self.assertRaises(WorkerFailed):await w.execute(b'',cmd)
  self.assertFalse(w.pids);self.assertEqual(w.admitted,0)
 async def test_unlimited_output_and_shutdown_reap(self):
  w=Workers();cmd=[sys.executable,'-c','import os\nwhile True: os.write(1,b"x"*65536)']
  with self.assertRaises(WorkerFailed):await asyncio.wait_for(w.execute(b'',cmd),5)
  self.assertFalse(w.pids)
  a=asyncio.create_task(w.execute(b'',command(10)));await asyncio.sleep(.03)
  b=asyncio.create_task(w.execute(b'',command(10)));await asyncio.sleep(.01)
  await w.close();await asyncio.gather(a,b,return_exceptions=True)
  self.assertFalse(w.pids);self.assertEqual(w.admitted,0);self.assertFalse(w.tasks)
  with self.assertRaises(Busy):await w.execute(b'',command())
 async def test_distinct_task_and_waiter_caps(self):
  with tempfile.TemporaryDirectory(dir=SCRATCH) as root:
   shared=Shared(DiskCache(root));gate=asyncio.Event()
   tasks=[asyncio.create_task(shared.run(str(i),gate.wait)) for i in range(3)]
   await asyncio.sleep(.01)
   with self.assertRaises(Busy):await shared.run('fourth',gate.wait)
   shared.max_waiters=2
   duplicate=asyncio.create_task(shared.run('0',gate.wait));await asyncio.sleep(.01)
   with self.assertRaises(Busy):await shared.run('0',gate.wait)
   gate.set();await asyncio.gather(*tasks,duplicate);self.assertFalse(shared.entries)
 async def test_deadline_kills_child_after_parent_exits(self):
  with tempfile.TemporaryDirectory(dir=SCRATCH) as root:
   record=Path(root)/'child.pid'
   script='import subprocess,sys; from pathlib import Path; p=subprocess.Popen([sys.executable,"-c","import time; time.sleep(20)"]); Path('+repr(str(record))+').write_text(str(p.pid)); print("parent exited")'
   w=Workers(deadline=.25)
   try:
    with self.assertRaises(Deadline):await asyncio.wait_for(w.execute(b'',[sys.executable,'-c',script]),3)
    pid=int(record.read_text());stat=Path(f'/proc/{pid}/stat')
    self.assertTrue(not stat.exists() or stat.read_text().split()[2]=='Z')
    self.assertFalse(w.pids);self.assertEqual(w.admitted,0)
   finally:
    if record.exists():
     try:os.kill(int(record.read_text()),9)
     except ProcessLookupError:pass
 async def test_shared_one_and_all_waiter_cancellation(self):
  with tempfile.TemporaryDirectory(dir=SCRATCH) as root:
   cache=DiskCache(root,diagnostics=True);shared=Shared(cache);w=Workers(diagnostics=True)
   factory=lambda:w.execute(b'',command(.15))
   a=asyncio.create_task(shared.run('same',factory));b=asyncio.create_task(shared.run('same',factory));await asyncio.sleep(.04)
   a.cancel()
   with self.assertRaises(asyncio.CancelledError):await a
   self.assertEqual(await b,b'ok');self.assertEqual(w.counters['launches'],1);self.assertFalse(shared.entries)
   a=asyncio.create_task(shared.run('all',lambda:w.execute(b'',command(10))));await asyncio.sleep(.04);a.cancel()
   with self.assertRaises(asyncio.CancelledError):await a
   self.assertFalse(w.pids);self.assertEqual(w.admitted,0);self.assertFalse(shared.entries)
   self.assertEqual(await shared.run('all',lambda:w.execute(b'',command(.01))),b'ok')
 async def test_disconnect_does_not_kill_other_waiter(self):
  class Client:
   def __init__(self):self.gone=False
   async def is_disconnected(self):return self.gone
  with tempfile.TemporaryDirectory(dir=SCRATCH) as root:
   cache=DiskCache(root);shared=Shared(cache);w=Workers(diagnostics=True);client=Client()
   factory=lambda:w.execute(b'',command(.2))
   a=asyncio.create_task(shared.run('key',factory,client));b=asyncio.create_task(shared.run('key',factory));await asyncio.sleep(.04);client.gone=True
   with self.assertRaises(asyncio.CancelledError):await a
   self.assertEqual(await b,b'ok');self.assertEqual(w.counters['launches'],1)


class CacheTests(unittest.TestCase):
 def test_entry_count_cap(self):
  with tempfile.TemporaryDirectory(dir=SCRATCH) as root:
   c=DiskCache(root,limit=1000);c.max_units=3
   for i in range(10):self.assertTrue(c.put(c.root/f'{i}.json',b'x'))
   self.assertEqual(len(c.units()),3)
 def test_canonical_all_query_and_body_fields(self):
  a=[('state','raw'),('ligand','A'),('future','1')]
  self.assertEqual(canonical_key('s','ligand',a,None),canonical_key('s','ligand',list(reversed(a)),None))
  self.assertNotEqual(canonical_key('s','ligand',a,None),canonical_key('s','ligand',a[:-1],None))
  self.assertNotEqual(canonical_key('s','ligand',[('a','1'),('a','2')],None),canonical_key('s','ligand',[('a','2'),('a','1')],None))
  self.assertEqual(canonical_key('s','analyze',[],{'x':1,'y':2}),canonical_key('s','analyze',[],{'y':2,'x':1}))
 def test_eviction_incoming_pins_and_atomic_analysis_pair(self):
  with tempfile.TemporaryDirectory(dir=SCRATCH) as root:
   c=DiskCache(root,limit=100,diagnostics=True);a=c.root/'a.json';b=c.root/'b.json'
   self.assertTrue(c.put(a,b'a'*60));c.protected.add(a)
   self.assertFalse(c.put(b,b'b'*50));self.assertEqual(c.read(a),b'a'*60)
   c.protected.clear();self.assertTrue(c.put(b,b'b'*50));self.assertFalse(a.exists());self.assertLessEqual(c.size(),100)
   with c.stage() as (folder,existing):
    (folder/'key.json').write_bytes(b'j'*20);(folder/'key.npy').write_bytes(b'n'*40);c.promote(folder,existing)
   self.assertTrue((c.analysis/'key.json').exists());self.assertTrue((c.analysis/'key.npy').exists());self.assertLessEqual(c.size(),100)
   self.assertTrue(c.put(a,b'a'*90));self.assertFalse((c.analysis/'key.json').exists());self.assertFalse((c.analysis/'key.npy').exists())
   self.assertFalse(list(c.root.rglob('.incoming-*')))
   print('EVICTION',json.dumps({'cache_bytes':c.size(),'limit':c.limit,'counters':c.counters}))
 def test_stage_keeps_in_use_evicted_bytes_and_no_partial_publish(self):
  with tempfile.TemporaryDirectory(dir=SCRATCH) as root:
   c=DiskCache(root,limit=100);target=c.analysis/'abc.npy';c.put(target,b'x'*60)
   with c.stage() as (folder,existing):
    self.assertTrue(c.put(c.root/'new.json',b'y'*70));self.assertFalse(target.exists())
    self.assertEqual((folder/'abc.npy').read_bytes(),b'x'*60)
    (folder/'partial.tmp').write_bytes(b'z'*40)
   self.assertFalse((c.analysis/'partial.tmp').exists());self.assertLessEqual(c.size(),100)


class ASGITests(unittest.IsolatedAsyncioTestCase):
 @classmethod
 def setUpClass(cls):
  cls.tmp=tempfile.TemporaryDirectory(prefix='runtime-lane-',dir=SCRATCH);cls.root=Path(cls.tmp.name)
  cls.modules=[]
  cls.source_hashes={str(p.relative_to(APP)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [APP/'manifest.json',APP/'backend/science.py',APP/'backend/main.py',APP/'backend/frozen.py',APP/'backend/runtime_worker.py']}
  cls.first=cls.copy('one');cls.second=cls.copy('two',variant=True)
 @classmethod
 def copy(cls,name,variant=False):
  root=cls.root/name;root.mkdir()
  shutil.copytree(APP/'backend',root/'backend',ignore=shutil.ignore_patterns('__pycache__'))
  shutil.copytree(APP/'static',root/'static')
  for file in ['manifest.json','structures.json']:shutil.copyfile(APP/file,root/file)
  # Copy reference packets through real paths; keep immutable science shared.
  (root/'data').mkdir()
  for source in (APP/'data').iterdir():
   if source.name=='revealed-reference':shutil.copytree(source,root/'data'/source.name)
   else:(root/'data'/source.name).symlink_to(source,target_is_directory=source.is_dir())
  (root/'pocket-residue-sets-20261028').symlink_to(APP/'pocket-residue-sets-20261028',target_is_directory=True)
  if variant:
   # Only isolated copies change, never active scientific source.
   with (root/'static/app.js').open('a') as f:f.write('\n// isolated revision test\n')
   with (root/'backend/science.py').open('a') as f:f.write('\n# isolated revision test\n')
  package=types.ModuleType('runtime_lane_'+name);package.__path__=[str(root/'backend')];sys.modules[package.__name__]=package
  with patch.dict(os.environ,{'CYP_DIAGNOSTICS':'1','CYP_RECEPTOR31_DIR':''}):module=importlib.import_module(package.__name__+'.runtime')
  cls.modules.append(module);return module
 @classmethod
 def tearDownClass(cls):
  for path,digest in cls.source_hashes.items():
   if hashlib.sha256((APP/path).read_bytes()).hexdigest()!=digest:
    print('CONCURRENT_SOURCE_CHANGE',path)
  cls.tmp.cleanup()
 async def asyncSetUp(self):
  asyncio.get_running_loop().set_debug(False)
  self.r=self.first
  self.client=httpx.AsyncClient(transport=httpx.ASGITransport(app=self.r.app),base_url='http://isolated')
 async def asyncTearDown(self):await self.client.aclose()
 async def test_headers_etags_static_and_error_allowlist(self):
  r=self.r
  boot=await self.client.get(r.PREFIX+'api/bootstrap');self.assertEqual(boot.headers['cache-control'],'no-store');self.assertNotIn('etag',boot.headers)
  b=boot.json();self.assertEqual(b['source_snapshot'],json.loads((APP/'manifest.json').read_text())['snapshot']);self.assertEqual(b['analysis_kernel_sha256'],r.REVISION['backend/science.py'])
  self.assertEqual(b['asset_policy'],'revalidate');self.assertEqual(b['asset_version'],r.SID)
  self.assertEqual(b['projection_version'],'pcoa-canonical-eigenspaces-v2');self.assertEqual(b['cluster_algorithm'],'average-linkage-original-distances-v1')
  for path in ['', 'static/app.js','static/style.css','static/residue-hbonds/index.json','static/vendor/molstar.js','static/vendor/molstar.css']:
   first=await self.client.get(r.PREFIX+path)
   self.assertEqual(first.status_code,200);self.assertEqual(first.headers['cache-control'],'no-cache')
   self.assertEqual(first.headers['etag'],'"'+hashlib.sha256(first.content).hexdigest()+'"')
   match=await self.client.get(r.PREFIX+path,headers={'If-None-Match':'"other", W/'+first.headers['etag']})
   self.assertEqual(match.status_code,304);self.assertEqual(match.content,b'')
  for path in ['/health','/unknown','/distances/'+'0'*24]:
   x=await self.client.get(r.BASE+path,headers={'If-None-Match':'*'});self.assertEqual(x.headers['cache-control'],'no-store');self.assertNotEqual(x.status_code,304)
  meta=await self.client.get(r.BASE+'/meta');self.assertIn('immutable',meta.headers['cache-control'])
  self.assertEqual((await self.client.get(r.BASE+'/meta',headers={'If-None-Match':meta.headers['etag']})).status_code,304)
  self.assertEqual((await self.client.post(r.BASE+'/meta')).status_code,405)
  unsupported=await self.client.request('CONNECT',r.BASE+'/meta')
  self.assertEqual(unsupported.status_code,405);self.assertEqual(unsupported.headers['cache-control'],'no-store')
 async def test_mixed_nested_snapshot_rejected(self):
  r=self.r
  for field in ('columns','analysis_spec','revealed_reference'):
   with self.assertRaises(r.WorkerFailed):r.bind_revision({'snapshot':'one',field:{'snapshot':'two'}})
  valid=r.bind_revision({'snapshot':'one','analysis_spec':{'snapshot':'one','source_snapshot':'science'},'columns':{'snapshot':'one'}})
  self.assertEqual(valid['snapshot'],r.SID);self.assertEqual(valid['columns']['snapshot'],r.SID);self.assertEqual(valid['analysis_spec']['snapshot'],r.SID)
  self.assertEqual(valid['analysis_spec']['source_snapshot'],'science');self.assertEqual(valid['columns']['source_snapshot'],'one')
 async def test_reference_files_are_revision_bound_and_changes_reject(self):
  r=self.copy('reference_mutation')
  path=r.ROOT/'data/revealed-reference/OCNT-2312554.json'
  self.assertIn('data/revealed-reference/manifest.json',r.REVISION)
  self.assertIn(str(path.relative_to(r.ROOT)),r.REVISION)
  self.assertTrue(r.active_revision())
  path.write_bytes(path.read_bytes()+b' ')
  self.assertFalse(r.active_revision())
  async with httpx.AsyncClient(transport=httpx.ASGITransport(app=r.app),base_url='http://isolated') as client:
   response=await client.get(r.PREFIX+'api/bootstrap')
   self.assertEqual(response.status_code,503);self.assertEqual(response.headers['cache-control'],'no-store')
 async def test_two_revisions_and_old_tabs(self):
  a=self.first;b=self.second
  self.assertNotEqual(a.SID,b.SID);self.assertEqual(a.SOURCE_SNAPSHOT,b.SOURCE_SNAPSHOT)
  async with httpx.AsyncClient(transport=httpx.ASGITransport(app=b.app),base_url='http://isolated') as client:
   old=await client.get(a.BASE+'/meta');self.assertEqual(old.status_code,404);self.assertIn('reload',old.json()['detail']);self.assertEqual(old.headers['cache-control'],'no-store')
   fresh=(await client.get(b.PREFIX+'api/bootstrap')).json();self.assertEqual(fresh['snapshot'],b.SID)
   old_asset=await self.client.get(a.PREFIX+'static/app.js')
   new_asset=await client.get(b.PREFIX+'static/app.js',headers={'If-None-Match':old_asset.headers['etag']})
   self.assertEqual(new_asset.status_code,200);self.assertNotEqual(old_asset.headers['etag'],new_asset.headers['etag'])
   self.assertNotEqual(a.REVISION['backend/science.py'],fresh['analysis_kernel_sha256'])
  print('REVISIONS',json.dumps({'one':a.SID,'two':b.SID,'source':a.SOURCE_SNAPSHOT}))
 async def test_real_worker_canonical_dedupe_and_analysis(self):
  r=self.r;before=r.workers.counters['launches'];url=r.BASE+'/ligand?ligand=OCNT-2395460&state=relaxed&lane=dedupe'
  start=time.monotonic()
  with patch.dict(os.environ,{'PYTHONPATH':str(r.ROOT)}):
   a,b=await asyncio.gather(self.client.get(url),self.client.get(r.BASE+'/ligand?lane=dedupe&state=relaxed&ligand=OCNT-2395460'))
  self.assertEqual(a.status_code,200,a.text[:200]);self.assertEqual(a.content,b.content);self.assertEqual(r.workers.counters['launches']-before,1)
  data=a.json();self.assertEqual(data['snapshot'],r.SID)
  self.assertEqual(data['analysis_spec']['snapshot'],r.SID);self.assertIn('source_snapshot',data['analysis_spec'])
  hit=await self.client.get(url);self.assertEqual(hit.headers['x-origin-cache'],'hit');self.assertEqual(hit.content,a.content)
  conditional=await self.client.get(url,headers={'If-None-Match':a.headers['etag']});self.assertEqual(conditional.status_code,304)
  pose_ids=[x['id'] for x in data['poses'][:3]]
  with patch.dict(os.environ,{'PYTHONPATH':str(r.ROOT)}):
   for state in ['raw','relaxed','paired']:
    ids=pose_ids if state!='paired' else [pid+'::raw' for pid in pose_ids]+[pid+'::relaxed' for pid in pose_ids]
    out=await self.client.post(r.BASE+'/analyze',json={'ligand':'OCNT-2395460','state':state,'metric':'rmsd','threshold':2,'atoms':[],'pose_ids':ids},headers={'If-None-Match':'*'})
    self.assertEqual(out.status_code,200,out.text[:200]);self.assertEqual(out.headers['cache-control'],'no-store');self.assertNotIn('etag',out.headers)
    analysis=out.json();self.assertEqual(analysis['snapshot'],r.SID)
    self.assertEqual(analysis['analysis_spec']['snapshot'],r.SID)
    matrix=await self.client.get(analysis['distances_url']);self.assertEqual(matrix.status_code,200)
    self.assertEqual((await self.client.get(analysis['distances_url'],headers={'If-None-Match':matrix.headers['etag']})).status_code,304)
  self.assertFalse(r.shared.entries);self.assertEqual(r.workers.admitted,0);self.assertLessEqual(r.cache.size(),r.MAX_CACHE)
  with patch.dict(os.environ,{'PYTHONPATH':str(r.ROOT)}):
   reused=await self.client.post(r.BASE+'/analyze',json={'ligand':'OCNT-2395460','state':'paired','metric':'rmsd','threshold':3,'atoms':[],'pose_ids':ids})
  self.assertEqual(reused.status_code,200,reused.text[:200]);self.assertTrue(reused.json()['diagnostics']['phase_cache_hit'])
  self.assertEqual(reused.json()['diagnostics']['distance_kernel_calls'],0)
  retry_body={'ligand':'OCNT-2395460','state':'paired','metric':'rmsd','threshold':3,'atoms':[],'pose_ids':ids}
  response_path=r.CACHE/(canonical_key(r.SID,'analyze',[],retry_body)+'.json');self.assertTrue(response_path.exists())
  distance_path=r.cache.analysis/(reused.json()['distances_url'].rsplit('/',1)[-1]+'.npy')
  limit=r.cache.limit;r.cache.protected.add(response_path)
  try:
   r.cache.limit=len(reused.content);r.cache.trim()
  finally:r.cache.limit=limit;r.cache.protected.discard(response_path)
  self.assertTrue(response_path.exists());self.assertFalse(distance_path.exists())
  launches=r.workers.counters['launches']
  with patch.dict(os.environ,{'PYTHONPATH':str(r.ROOT)}):retried=await self.client.post(r.BASE+'/analyze',json=retry_body)
  self.assertEqual(retried.status_code,200,retried.text[:200]);self.assertEqual(r.workers.counters['launches']-launches,1)
  self.assertEqual((await self.client.get(retried.json()['distances_url'])).status_code,200)
  print('REAL_WORKER',json.dumps({'seconds':time.monotonic()-start,'ligand_bytes':len(a.content),'worker_counters':r.workers.counters,'cache_counters':r.cache.counters,'cache_bytes':r.cache.size(),'origin':a.headers.get('x-origin-cache'),'timing':a.headers.get('server-timing')}))
 async def test_real_background_foreground_progress_and_distinct_keys(self):
  r=self.r;before=r.workers.counters['launches'];order=[]
  async def request(path,label):
   response=await self.client.get(r.BASE+path);order.append(label);return response
  with patch.dict(os.environ,{'PYTHONPATH':str(r.ROOT)}):
   bg=asyncio.create_task(request('/receptor-background?ligand=OCNT-2313300&state=relaxed&lane=progress','background'))
   await asyncio.sleep(.05)
   rejected=await self.client.get(r.BASE+'/receptor-background?ligand=OCNT-2312554&state=raw&lane=progress')
   self.assertEqual(rejected.status_code,429);self.assertEqual(rejected.headers['cache-control'],'no-store')
   foreground=asyncio.create_task(request('/ligand?ligand=OCNT-2312554&state=relaxed&format=columns-v1&lane=progress','foreground'))
   a,b=await asyncio.gather(bg,foreground)
  self.assertEqual(a.status_code,200,a.text[:200]);self.assertEqual(b.status_code,200,b.text[:200])
  self.assertEqual(r.workers.counters['launches']-before,2);self.assertFalse(r.shared.entries)
  print('REAL_PROGRESS',json.dumps({'order':order,'launches':r.workers.counters['launches']-before,'background_bytes':len(a.content),'foreground_bytes':len(b.content)}))
 async def test_real_eviction_and_retry(self):
  r=self.r;url=r.BASE+'/receptor-background?ligand=OCNT-2313300&state=relaxed&lane=eviction';before=r.workers.counters['launches']
  with patch.dict(os.environ,{'PYTHONPATH':str(r.ROOT)}):
   first=await self.client.get(url);self.assertEqual(first.status_code,200)
   limit=r.cache.limit;evictions=r.cache.counters['evictions']
   try:
    r.cache.limit=1;r.cache.trim();self.assertLessEqual(r.cache.size(),1)
   finally:r.cache.limit=limit
   second=await self.client.get(url);self.assertEqual(second.status_code,200)
  self.assertEqual(first.content,second.content);self.assertEqual(r.workers.counters['launches']-before,2)
  self.assertGreater(r.cache.counters['evictions'],evictions);self.assertLessEqual(r.cache.size(),r.MAX_CACHE)
  print('REAL_EVICTION_RETRY',json.dumps({'launches':r.workers.counters['launches']-before,'evictions':r.cache.counters['evictions']-evictions,'restored_limit':r.cache.limit,'cache_bytes':r.cache.size()}))
 async def test_real_paired_compact_worker_limit(self):
  r=self.r
  cases=[('OCNT-2312554','paired'),('OCNT-2313300','paired')]
  for ligand,state in cases:
   with patch.dict(os.environ,{'PYTHONPATH':str(r.ROOT)}):
    response=await self.client.get(r.BASE+f'/ligand?ligand={ligand}&state={state}&format=columns-v1&lane=largest')
   self.assertEqual(response.status_code,200,response.text[:200]);data=response.json()
   self.assertEqual(data['snapshot'],r.SID);self.assertEqual(data['analysis_spec']['snapshot'],r.SID)
   self.assertEqual(data['columns']['snapshot'],r.SID)
   self.assertLess(len(response.content),32*1024**2);self.assertLessEqual(r.cache.size(),r.MAX_CACHE)
   print('REAL_PAIRED',json.dumps({'ligand':ligand,'state':state,'bytes':len(response.content),'cache_bytes':r.cache.size(),'snapshot':data['snapshot']}))
 async def test_coordinate_gzip_content_identity_and_304(self):
  r=self.r;key,entry=next(iter(r.structures.items()))
  root=Path(os.environ['CYP_COORDINATE_ROOT']);compressed=(root/entry['path']).read_bytes()
  self.assertEqual(hashlib.sha256(gzip.decompress(compressed)).hexdigest(),entry['sha256'])
  with patch.object(r,'COORDS',root):
   result=await self.client.get(r.PREFIX+entry['path']);self.assertEqual(result.status_code,200)
   self.assertEqual(result.content,compressed);self.assertEqual(result.headers['cache-control'],'no-cache');self.assertNotIn('content-encoding',result.headers)
   self.assertEqual(result.headers['etag'],'"'+hashlib.sha256(compressed).hexdigest()+'"')
   conditional=await self.client.get(r.PREFIX+entry['path'],headers={'If-None-Match':result.headers['etag']});self.assertEqual(conditional.status_code,304)
  print('COORDINATE',json.dumps({'key':key,'declared_sha256':entry['sha256'],'gzip_sha256':hashlib.sha256(compressed).hexdigest(),'gzip_bytes':len(compressed)}))
 async def test_default_diagnostics_disabled(self):
  root=self.root/'one';package=types.ModuleType('runtime_lane_default');package.__path__=[str(root/'backend')];sys.modules[package.__name__]=package
  with patch.dict(os.environ,{'CYP_DIAGNOSTICS':'0','CYP_RECEPTOR31_DIR':''}):module=importlib.import_module(package.__name__+'.runtime')
  self.assertFalse(module.DIAGNOSTICS)
  async with httpx.AsyncClient(transport=httpx.ASGITransport(app=module.app),base_url='http://isolated') as client:
   response=await client.get(module.BASE+'/diagnostics');self.assertEqual(response.status_code,404);self.assertEqual(response.headers['cache-control'],'no-store')
 async def test_revision_change_detected_without_serving_stale_immutable(self):
  r=self.second;path=r.ROOT/'static/style.css';original=path.read_bytes()
  try:
   path.write_bytes(original+b'\n/* isolated write */\n')
   async with httpx.AsyncClient(transport=httpx.ASGITransport(app=r.app),base_url='http://isolated') as client:
    for target in [r.PREFIX+'api/bootstrap',r.BASE+'/meta',r.PREFIX+'static/app.js']:
     result=await client.get(target);self.assertEqual(result.status_code,503);self.assertEqual(result.headers['cache-control'],'no-store')
  finally:
   path.write_bytes(original)
   # Restore the private test copy stat receipt for later tests only.
   r.REVISION_STATS['static/style.css']=(path.stat().st_mtime_ns,path.stat().st_size,path.stat().st_ino)
 async def test_http_failure_timeout_overload_retry(self):
  r=self.r;execute=r.workers.execute;mode='fail'
  async def injected(payload,**kwargs):
   if mode=='fail':return await execute(payload,command=command(.01,True))
   if mode=='timeout':raise r.Deadline()
   if mode=='busy':raise r.Busy()
   return json.dumps({'status':200,'result':{'fixture':True,'snapshot':'source','analysis_spec':{'snapshot':'source'}}}).encode()
  with patch.object(r.workers,'execute',injected):
   for mode,status in [('fail',503),('timeout',504),('busy',429),('ok',200)]:
    x=await self.client.get(r.BASE+'/ligand?ligand=fixture&test=failure')
    self.assertEqual(x.status_code,status,x.text)
    if status!=200:self.assertEqual(x.headers['cache-control'],'no-store');self.assertNotIn('etag',x.headers)
    if status==429:self.assertEqual(x.headers['retry-after'],'1')
    self.assertFalse(r.shared.entries);self.assertEqual(r.workers.admitted,0)
  self.assertEqual(x.json()['analysis_spec']['snapshot'],r.SID)
 async def test_asgi_waiter_cancellation_and_retry(self):
  r=self.r;execute=r.workers.execute;before=r.workers.counters['launches']
  fixture=json.dumps({'status':200,'result':{'snapshot':'fixture','fixture':True}})
  async def injected(payload,**kwargs):return await execute(payload,command=command(.25,value=fixture))
  with patch.object(r.workers,'execute',injected):
   url=r.BASE+'/ligand?ligand=fixture&lane=one-cancel'
   a=asyncio.create_task(self.client.get(url));b=asyncio.create_task(self.client.get(url));await asyncio.sleep(.06)
   a.cancel()
   with self.assertRaises(asyncio.CancelledError):await a
   response=await b;self.assertEqual(response.status_code,200);self.assertEqual(r.workers.counters['launches']-before,1)
   url=r.BASE+'/ligand?ligand=fixture&lane=all-cancel'
   a=asyncio.create_task(self.client.get(url));b=asyncio.create_task(self.client.get(url));await asyncio.sleep(.06)
   a.cancel();b.cancel();await asyncio.gather(a,b,return_exceptions=True)
   self.assertFalse(r.shared.entries);self.assertFalse(r.workers.pids);self.assertEqual(r.workers.admitted,0)
   retry=await self.client.get(url);self.assertEqual(retry.status_code,200)
   self.assertEqual(r.workers.counters['launches']-before,3)

if __name__=='__main__':unittest.main(verbosity=2)
