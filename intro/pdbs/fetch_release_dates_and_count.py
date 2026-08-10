#!/usr/bin/env python3
from pathlib import Path
import csv, json, requests
from datetime import date

root=Path('/home/dan/projects/ADMET-CYP/intro/pdbs')
rows=list(csv.DictReader((root/'manifest.csv').open()))
ids=sorted({r['pdb_id'].upper() for r in rows})
q='''query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id rcsb_accession_info{deposit_date initial_release_date}}}'''
r=requests.post('https://data.rcsb.org/graphql',json={'query':q,'variables':{'ids':ids}},timeout=120)
r.raise_for_status(); payload=r.json()
if payload.get('errors'): raise RuntimeError(payload['errors'])
entries=payload['data']['entries']; by={e['rcsb_id'].upper():e['rcsb_accession_info'] for e in entries}
cutoffs={'chai1':'2021-01-12','boltz2':'2023-06-01','protenix':'2021-09-30','openfold3_preview':'2021-09-30'}
out=[]
for x in rows:
 info=by.get(x['pdb_id'].upper())
 if not info: raise RuntimeError('missing RCSB date '+x['pdb_id'])
 y=dict(x); y.update(info)
 for model,c in cutoffs.items(): y['post_'+model]=str(info['initial_release_date'][:10]>c).lower()
 out.append(y)
with (root/'manifest_with_release_dates.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys()); w.writeheader(); w.writerows(out)
summary={}
for scope,subset in [('all_panel',out),('CYP3A4',[x for x in out if x['gene']=='CYP3A4'])]:
 methods={}
 for model,c in cutoffs.items():
  post=[x for x in subset if x['initial_release_date'][:10]>c]
  methods[model]={'cutoff':c,'post_cutoff_count':len(post),'total':len(subset),'pdb_ids':[x['pdb_id'] for x in post]}
 summary[scope]={'total':len(subset),'methods':methods}
(root/'post_training_cutoff_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
