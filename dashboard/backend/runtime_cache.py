"""Bounded serialized-response/disk-analysis cache and shared request tasks.

The retained cache is disk-only (128 MiB); no parsed response graph is retained.
Eviction is least-recently-used by file mtime. Analysis JSON/NPY pairs are one
unit. A read returns owned bytes before eviction can run. Temporary worker
outputs are isolated from published cache entries and promoted with reservation.
"""
import asyncio
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time
from contextlib import contextmanager
from .workers import Busy


def canonical_key(snapshot,operation,query,body):
 # Query names are sorted, but repeated values keep their order. All parameters
 # remain part of identity, including currently unknown ones.
 pairs=sorted(query,key=lambda item:item[0])
 value=[snapshot,operation,pairs,body]
 return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


class DiskCache:
 def __init__(self,root,limit=128*1024**2,diagnostics=False):
  self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
  self.analysis=self.root/'analysis';self.analysis.mkdir(exist_ok=True)
  self.limit=limit;self.max_units=512;self.diagnostics=diagnostics;self.protected=set()
  self.counters=dict(hits=0,misses=0,evictions=0,serialization_seconds=0.,shared_waiters=0)
  self.trim()
 def count(self,key,value=1):
  if self.diagnostics:self.counters[key]+=value
 def units(self):
  units={}
  for path in self.root.rglob('*'):
   if path.is_file():
    # Do not evict one half of a published scientific analysis pair.
    key=path.with_suffix('') if path.parent==self.analysis else path
    units.setdefault(key,[]).append(path)
  return units
 def size(self):return sum(p.stat().st_size for ps in self.units().values() for p in ps)
 def trim(self,incoming=0,exclude=(),incoming_units=0):
  units=self.units();size=sum(p.stat().st_size for ps in units.values() for p in ps)
  count=len(units)
  for key,paths in sorted(units.items(),key=lambda item:min(p.stat().st_mtime_ns for p in item[1])):
   if size+incoming<=self.limit and count+incoming_units<=self.max_units:break
   if any(p in self.protected or p in exclude for p in paths):continue
   for p in paths:size-=p.stat().st_size;p.unlink(missing_ok=True)
   count-=1
   self.count('evictions')
  return size+incoming<=self.limit and count+incoming_units<=self.max_units
 def read(self,path):
  try:
   payload=path.read_bytes();os.utime(path,None);self.count('hits');return payload
  except FileNotFoundError:self.count('misses');return None
 def put(self,path,payload):
  if len(payload)>self.limit:return False
  if path.exists():return True
  unit=path.with_suffix('') if path.parent==self.analysis else path
  incoming_units=int(unit not in self.units())
  if not self.trim(len(payload),incoming_units=incoming_units):return False
  path.parent.mkdir(parents=True,exist_ok=True)
  fd,name=tempfile.mkstemp(prefix='.incoming-',dir=path.parent)
  try:
   with os.fdopen(fd,'wb') as out:out.write(payload);out.flush();os.fsync(out.fileno())
   os.replace(name,path)
  finally:Path(name).unlink(missing_ok=True)
  return True
 @contextmanager
 def stage(self,seed=True):
  # Hard links retain already-cached scientific kernels without copying bytes.
  # Newly computed files stay private until the worker exits successfully.
  scratch=Path(os.environ.get('TMPDIR',str(self.root.parent/'cache-work')))
  scratch.mkdir(parents=True,exist_ok=True)
  with tempfile.TemporaryDirectory(prefix='cyp-analysis-',dir=scratch) as folder:
   folder=Path(folder)
   existing={}
   if seed:
    for p in self.analysis.iterdir():
     if p.is_file() and p.suffix in {'.json','.npy','.npz'}:
      os.link(p,folder/p.name);existing[p.name]=(p.stat().st_ino,p.stat().st_size)
   yield folder,existing
 def promote(self,folder,existing):
  groups={}
  for p in folder.iterdir():
   if p.suffix not in {'.json','.npy','.npz'}:continue
   if existing.get(p.name)==(p.stat().st_ino,p.stat().st_size):
    # The worker has finished. Drop unused seed links before any eviction, so
    # evicted cache bytes do not remain retained outside the cache directory.
    p.unlink();continue
   groups.setdefault(p.stem,[]).append(p)
  for paths in groups.values():
   incoming=sum(p.stat().st_size for p in paths)
   targets=[self.analysis/p.name for p in paths]
   incoming_units=int(targets[0].with_suffix('') not in self.units())
   if incoming>self.limit or not self.trim(incoming,exclude=targets,incoming_units=incoming_units):continue
   # Both names publish without an await; an HTTP reader cannot interleave.
   for p,target in zip(paths,targets):
    if target.exists():continue
    self.put(target,p.read_bytes())


class Shared:
 def __init__(self,cache):self.entries={};self.cache=cache;self.max_entries=3;self.max_waiters=64
 async def run(self,key,factory,request=None):
  entry=self.entries.get(key)
  if entry is None:
   if len(self.entries)>=self.max_entries:raise Busy('Distinct request limit; retry later')
   entry={'task':asyncio.create_task(factory()),'waiters':0};self.entries[key]=entry
  else:
   if entry['waiters']>=self.max_waiters:raise Busy('Shared waiter limit; retry later')
   self.cache.count('shared_waiters')
  entry['waiters']+=1
  task=entry['task'];watch=None
  async def disconnected():
   assert request is not None
   while True:
    if await request.is_disconnected():return
    await asyncio.sleep(.05)
  try:
   if request is None:return await asyncio.shield(task)
   watch=asyncio.create_task(disconnected())
   done,_=await asyncio.wait((task,watch),return_when=asyncio.FIRST_COMPLETED)
   if task in done:return task.result()
   raise asyncio.CancelledError('Request disconnected')
  finally:
   if watch is not None:
    watch.cancel();await asyncio.gather(watch,return_exceptions=True)
   entry['waiters']-=1
   if not entry['waiters']:
    # Remove before awaiting cleanup so a retry cannot join a cancelled job.
    if self.entries.get(key) is entry:self.entries.pop(key)
    if not task.done():task.cancel()
    await asyncio.gather(task,return_exceptions=True)


def worker_main():
 """Use the active scientific worker unchanged, with private cache writes."""
 import runpy
 from .frozen import FrozenStore
 folder=Path(os.environ['CYP_RUNTIME_ANALYSIS_STAGE'])
 analysis=Path(os.environ['CYP_RUNTIME_ANALYSIS_CACHE'])
 existing={}
 # Seed only after admission, inside the one active child. Queued jobs hold
 # empty working directories, not old cached analysis files.
 for p in analysis.iterdir():
  if p.is_file() and p.suffix in {'.json','.npy','.npz'}:
   os.link(p,folder/p.name);existing[p.name]=(p.stat().st_ino,p.stat().st_size)
 (folder/'.seeds').write_text(json.dumps(existing))
 original=FrozenStore.__init__
 def staged(self,*args,**kwargs):
  original(self,*args,**kwargs)
  self.analysis_folder=folder
 FrozenStore.__init__=staged
 runpy.run_module('backend.runtime_worker',run_name='__main__')

if __name__=='__main__':worker_main()
