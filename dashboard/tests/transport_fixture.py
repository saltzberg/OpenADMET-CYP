#!/usr/bin/env python3
"""Small transport test fixture from current authorized frozen source rows."""
import json,os,sys
from pathlib import Path
sys.path.insert(0,str(Path(os.environ.get('CYP_TEST_APP', Path(__file__).resolve().parents[1])).resolve()))
from backend.frozen import FrozenStore
from backend.compact import encode
import tempfile
with tempfile.TemporaryDirectory(prefix='transport-contract-') as cache:
 s=FrozenStore('OCNT-2312554');s.analysis_folder=Path(cache)
 data=s.ligand('OCNT-2312554','raw');data['poses']=data['poses'][:3]
# This is a codec fixture, not a changed published population or science result.
data['coverage']={'shown':3,'fixture':True}
print(json.dumps({'source':data,'wire':encode(data)},separators=(',',':'),allow_nan=False))
