"""Lossless columns-v1 ligand transport. No numeric conversion or rounding.

Only scores, percentiles and contacts move out of poses. Other fields stay exact.
Consumers can read columns synchronously without decoding the legacy row graph.
"""
import hashlib
import json

FORMAT='columns-v1'

def population_sha256(ids):
 return hashlib.sha256(json.dumps(ids,separators=(',',':')).encode()).hexdigest()

def encode(data):
 poses=data['poses'];ids=[p['id'] for p in poses]
 if len(set(ids))!=len(ids):raise ValueError('Duplicate pose identity')
 keys={field:list(dict.fromkeys(k for p in poses for k in p.get(field,{}))) for field in ('scores','percentiles')}
 vocab=list(dict.fromkeys(c for p in poses for c in p.get('contacts',[])));lookup={c:i for i,c in enumerate(vocab)}
 columns={field:[[p.get(field,{}).get(k) for p in poses] for k in keys[field]] for field in keys}
 # Presence masks retain absent maps/keys separately from explicit null values.
 presence={field:[None if all(k in p.get(field,{}) for p in poses) else [k in p.get(field,{}) for p in poses] for k in keys[field]] for field in keys}
 result={k:v for k,v in data.items() if k!='poses'}
 result.update(format=FORMAT,poses=[{k:v for k,v in p.items() if k not in ('scores','percentiles','contacts')} for p in poses],columns=dict(version=1,snapshot=data['snapshot'],ligand=data['ligand'],state=data['state'],pose_ids=ids,population_sha256=population_sha256(ids),dimensions=dict(poses=len(poses),scores=len(keys['scores']),percentiles=len(keys['percentiles']),contacts=len(vocab)),score_keys=keys['scores'],percentile_keys=keys['percentiles'],scores=columns['scores'],percentiles=columns['percentiles'],presence=presence,map_presence={f:[f in p for p in poses] for f in ('scores','percentiles','contacts')},contact_vocab=vocab,contacts=[[lookup[c] for c in p.get('contacts',[])] for p in poses]))
 return result

def validate(data):
 if data.get('format')!=FORMAT:raise ValueError('Unknown compact format')
 c=data.get('columns',{})
 if c.get('version')!=1:raise ValueError('Unknown columns version')
 for k in ('snapshot','ligand','state'):
  if c.get(k)!=data.get(k):raise ValueError('Mixed compact revision or population')
 ids=c.get('pose_ids');poses=data.get('poses')
 if isinstance(poses,list) and any(not isinstance(p,dict) or any(field in p for field in ('scores','percentiles','contacts')) for p in poses):raise ValueError('Mixed legacy fields in compact pose row')
 if not isinstance(ids,list) or not isinstance(poses,list) or ids!=[p.get('id') for p in poses] or len(set(ids))!=len(ids):raise ValueError('Invalid ordered pose identities')
 if population_sha256(ids)!=c.get('population_sha256'):raise ValueError('Population hash mismatch')
 n=len(ids);dims=c.get('dimensions',{})
 if dims.get('poses')!=n:raise ValueError('Invalid pose dimension')
 for field,key in (('scores','score_keys'),('percentiles','percentile_keys')):
  names=c.get(key,[]);values=c.get(field,[]);masks=c.get('presence',{}).get(field,[])
  if len(set(names))!=len(names) or dims.get(field)!=len(names) or len(values)!=len(names) or len(masks)!=len(names):raise ValueError('Invalid metric dimension')
  if any(not isinstance(a,list) or len(a)!=n for a in values) or any(a is not None and (not isinstance(a,list) or len(a)!=n) for a in masks):raise ValueError('Invalid metric row dimension')
  if any(type(v) is not bool for a in masks if a is not None for v in a):raise ValueError('Invalid metric presence')
 vocab=c.get('contact_vocab',[]);contacts=c.get('contacts',[])
 if len(set(vocab))!=len(vocab) or any(not isinstance(v,str) for v in vocab) or dims.get('contacts')!=len(vocab) or len(contacts)!=n:raise ValueError('Invalid contact dimension')
 if any(not isinstance(a,list) or any(type(i) is not int or i<0 or i>=len(vocab) for i in a) for a in contacts):raise ValueError('Invalid contact vocabulary index')
 for field in ('scores','percentiles','contacts'):
  mask=c.get('map_presence',{}).get(field,[])
  if len(mask)!=n or any(type(v) is not bool for v in mask):raise ValueError('Invalid map presence')
 return c

def decode(data):
 c=validate(data);out={k:v for k,v in data.items() if k not in ('format','columns','poses')};poses=[]
 for i,p in enumerate(data['poses']):
  p=dict(p)
  for field,key in (('scores','score_keys'),('percentiles','percentile_keys')):
   if c['map_presence'][field][i]:p[field]={k:c[field][j][i] for j,k in enumerate(c[key]) if c['presence'][field][j] is None or c['presence'][field][j][i]}
  if c['map_presence']['contacts'][i]:p['contacts']=[c['contact_vocab'][j] for j in c['contacts'][i]]
  poses.append(p)
 out['poses']=poses
 return out
