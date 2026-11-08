"""One process, two waiting jobs total, at most one background job.

Foreground jobs run first. The deadline includes queue time. Cancellation kills
and reaps the process group before releasing its slot. stdout/stderr are bounded
while reading, not after communicate has allocated an arbitrary response.
"""
import asyncio, os, signal, sys, time
from collections import deque

class Busy(Exception):pass
class Deadline(Exception):pass
class WorkerFailed(Exception):pass

class Workers:
 def __init__(self,concurrency=1,queued=2,deadline=45.,diagnostics=False):
  if concurrency!=1 or queued!=2:raise ValueError('Supported bounds: one active, two queued total')
  self.capacity=3;self.admitted=0;self.deadline=deadline;self.pids=set();self.last_killed=None
  self.active=None;self.foreground=deque();self.background=deque();self.background_admitted=0
  self.tasks=set();self.closed=False
  self.diagnostics=diagnostics
  self.counters=dict(launches=0,queue_seconds=0.,worker_seconds=0.,failures=0,timeouts=0,cancellations=0)
 def _count(self,key,value=1.):
  if self.diagnostics:self.counters[key]+=value
 def _dispatch(self):
  if self.active is not None:return
  for queue in (self.foreground,self.background):
   while queue:
    ticket=queue.popleft()
    if ticket.cancelled():continue
    self.active=ticket;ticket.set_result(None);return
 async def execute(self,payload,command=None,background=False,timing=None,env=None):
  # One background can never exhaust both waiting slots. Keep the original
  # total admission bound; priority does not require an extra queue slot.
  if self.closed or self.admitted>=self.capacity or (background and self.background_admitted):
   raise Busy('Queue full; retry later')
  ticket=asyncio.get_running_loop().create_future();queue=self.background if background else self.foreground
  queue.append(ticket);self.admitted+=1;self.background_admitted+=int(background)
  owner=asyncio.current_task();self.tasks.add(owner)
  process=None;start=time.monotonic();running=None
  try:
   self._dispatch()
   async with asyncio.timeout(self.deadline):
    await ticket
    running=time.monotonic();wait=running-start;self._count('queue_seconds',wait)
    if timing is not None:timing['queue']=wait
    async def launch():
     child=await asyncio.create_subprocess_exec(*(command or [sys.executable,'-m','backend.worker']),stdin=asyncio.subprocess.PIPE,stdout=asyncio.subprocess.PIPE,stderr=asyncio.subprocess.PIPE,start_new_session=True,env=dict(os.environ,**(env or {}),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1'))
     self.pids.add(child.pid);self._count('launches');return child
    spawning=asyncio.create_task(launch())
    try:process=await asyncio.shield(spawning)
    except asyncio.CancelledError:
     # Obtain the handle even if the client cancels during OS process creation.
     process=await spawning;raise
    async def bounded(stream,limit):
     out=bytearray()
     while True:
      chunk=await stream.read(65536)
      if not chunk:return bytes(out)
      out.extend(chunk)
      if len(out)>limit:raise WorkerFailed('Worker output exceeded limit')
    async def feed():
     assert process.stdin is not None
     try:process.stdin.write(payload);await process.stdin.drain()
     except (BrokenPipeError,ConnectionResetError):pass
     finally:process.stdin.close()
    tasks=[asyncio.create_task(feed()),asyncio.create_task(bounded(process.stdout,32*1024**2)),asyncio.create_task(bounded(process.stderr,65536))]
    try:
     _,stdout,_=await asyncio.gather(*tasks);await process.wait()
     if process.returncode:raise WorkerFailed('Scientific process exited unsuccessfully')
     return stdout
    finally:
     for task in tasks:
      if not task.done():task.cancel()
     await asyncio.gather(*tasks,return_exceptions=True)
  except TimeoutError:
   self._count('timeouts');raise Deadline('Analysis exceeded execution/queue deadline')
  except asyncio.CancelledError:self._count('cancellations');raise
  except WorkerFailed:self._count('failures');raise
  finally:
   if process:
    # A parent can exit while a child still owns its pipes. Kill the owned
    # process group even when the parent's returncode is already available.
    try:
     os.killpg(process.pid,signal.SIGKILL)
     self.last_killed=process.pid
    except ProcessLookupError:pass
    await process.communicate()
    self.pids.discard(process.pid)
   if running is not None:
    seconds=time.monotonic()-running;self._count('worker_seconds',seconds)
    if timing is not None:timing['worker']=seconds
   if self.active is ticket:self.active=None
   else:
    try:queue.remove(ticket)
    except ValueError:pass
   self.admitted-=1;self.background_admitted-=int(background);self._dispatch()
   self.tasks.discard(owner)
 async def close(self):
  self.closed=True
  tasks=list(self.tasks)
  for task in tasks:
   if task is not None:task.cancel()
  await asyncio.gather(*tasks,return_exceptions=True)
