"""Single-flight, last-good publication reader. Never opens producer files for writing."""
import hashlib,io,json,threading
import pandas as pd

class EnergySnapshot:
 def __init__(self,folder,validate):
  self.folder=folder;self.validate=validate;self.lock=threading.Lock()
  self.current=({},None);self.accepted=None;self.attempted=None;self.running=False;self.error=None;self.validations=0
 def generation(self):
  try:
   return tuple((s.st_ino,s.st_mtime_ns,s.st_ctime_ns,s.st_size) for s in [(self.folder/n).stat() for n in ('energies.parquet','energies.receipt.json')])
  except FileNotFoundError:return None
 def read(self):
  generation=self.generation()
  with self.lock:
   if generation!=self.attempted and not self.running:
    self.running=True;self.attempted=generation
    threading.Thread(target=self.update,args=(generation,),daemon=True,name='energy-snapshot-validation').start()
   values,version=self.current
   status='refreshing' if self.running else ('invalid' if self.error else ('verified' if version else 'unavailable'))
   return values,dict(version=version,status=status,rows=len(values),last_good=bool(version),error=self.error,validations=self.validations)
 def update(self,generation):
  try:
   if generation is None:raise ValueError('Complete energy table and receipt unavailable')
   receipt_bytes=(self.folder/'energies.receipt.json').read_bytes()
   payload=(self.folder/'energies.parquet').read_bytes()
   receipt=json.loads(receipt_bytes)
   digest=hashlib.sha256(payload).hexdigest()
   if receipt.get('sha256')!=digest:raise ValueError('Energy receipt hash mismatch')
   frame=pd.read_parquet(io.BytesIO(payload))
   if not frame.pose_id.is_unique or len(frame)!=receipt.get('rows'):raise ValueError('Energy receipt row count or identity mismatch')
   counts=frame.energy_status.value_counts().to_dict()
   if counts!=receipt.get('counts'):raise ValueError('Energy receipt status counts mismatch')
   values=self.validate(frame)
   if self.generation()!=generation or (self.folder/'energies.receipt.json').read_bytes()!=receipt_bytes:raise ValueError('Energy publication changed during validation')
   version=hashlib.sha256(receipt_bytes+payload).hexdigest()
   with self.lock:self.current=(values,version);self.accepted=generation;self.error=None
  except Exception as exc:
   with self.lock:self.error=f'{type(exc).__name__}: {exc}'
  finally:
   with self.lock:self.running=False;self.validations+=1
