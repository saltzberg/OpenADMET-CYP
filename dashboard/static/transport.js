/* Lossless columns-v1 accessors. Metric maps are views, not retained copies. */
(function(root,factory){'use strict';const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;else root.PoseTransport=api;})(globalThis,function(){'use strict';
const crypto={subtle:{digest:(algorithm,bytes)=>(globalThis.PoseHash||require('./sha256.js')).digest(algorithm,bytes)}};
const fail=message=>{throw new Error('Invalid compact ligand response: '+message);};
const same=(a,b)=>Array.isArray(a)&&Array.isArray(b)&&a.length===b.length&&a.every((x,i)=>x===b[i]);
async function validate(data){
 if(data?.format!=='columns-v1')fail('unknown format');
 const c=data.columns;if(c?.version!==1)fail('unknown schema version');
 for(const key of ['snapshot','ligand','state'])if(c[key]!==data[key])fail('mixed '+key);
 if(Array.isArray(data.poses)&&data.poses.some(p=>!p||typeof p!=='object'||Array.isArray(p)||['scores','percentiles','contacts'].some(field=>Object.hasOwn(p,field))))fail('mixed legacy fields in compact pose row');
 if(!Array.isArray(data.poses)||!Array.isArray(c.pose_ids)||!same(c.pose_ids,data.poses.map(p=>p.id))||new Set(c.pose_ids).size!==c.pose_ids.length)fail('ordered pose identity');
 const bytes=new TextEncoder().encode(JSON.stringify(c.pose_ids));

 const hash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes)),v=>v.toString(16).padStart(2,'0')).join('');
 if(hash!==c.population_sha256)fail('population hash');
 const n=c.pose_ids.length;if(c.dimensions?.poses!==n)fail('pose dimension');
 for(const[field,namesKey]of [['scores','score_keys'],['percentiles','percentile_keys']]){
  const keys=c[namesKey],values=c[field],masks=c.presence?.[field];
  if(!Array.isArray(keys)||keys.some(k=>typeof k!=='string')||new Set(keys).size!==keys.length||c.dimensions[field]!==keys.length||!Array.isArray(values)||values.length!==keys.length||!Array.isArray(masks)||masks.length!==keys.length)fail(field+' metric dimension');
  if(values.some(a=>!Array.isArray(a)||a.length!==n)||masks.some(a=>a!==null&&(!Array.isArray(a)||a.length!==n||a.some(v=>typeof v!=='boolean'))))fail(field+' row dimension/presence');
 }
 const vocab=c.contact_vocab,contacts=c.contacts;
 if(!Array.isArray(vocab)||vocab.some(v=>typeof v!=='string')||new Set(vocab).size!==vocab.length||c.dimensions.contacts!==vocab.length||!Array.isArray(contacts)||contacts.length!==n)fail('contact dimensions');
 if(contacts.some(a=>!Array.isArray(a)||a.some(i=>!Number.isInteger(i)||i<0||i>=vocab.length)))fail('contact vocabulary index');
 for(const field of ['scores','percentiles','contacts']){const mask=c.map_presence?.[field];if(!Array.isArray(mask)||mask.length!==n||mask.some(v=>typeof v!=='boolean'))fail(field+' map presence');}
 return c;
}
function mapView(c,field,keys,indexes,i){
 const present=j=>c.presence[field][j]===null||c.presence[field][j][i];
 return new Proxy(Object.create(null),{
  get(_target,key){const j=indexes.get(key);return j!==undefined&&present(j)?c[field][j][i]:undefined;},
  has(_target,key){const j=indexes.get(key);return j!==undefined&&present(j);},
  ownKeys(){return keys.filter((_key,j)=>present(j));},
  getOwnPropertyDescriptor(_target,key){const j=indexes.get(key);return j!==undefined&&present(j)?{enumerable:true,configurable:true,writable:false,value:c[field][j][i]}:undefined;},
  set(){throw new TypeError('Metric columns are read-only');},deleteProperty(){throw new TypeError('Metric columns are read-only');},defineProperty(){throw new TypeError('Metric columns are read-only');}
 });
}
async function adapt(data){
 // Legacy is supported for explicit compatibility; the application requests columns-v1.
 if(data?.format===undefined&&!data?.columns)return data;
 const c=await validate(data),indexes={scores:new Map(c.score_keys.map((k,i)=>[k,i])),percentiles:new Map(c.percentile_keys.map((k,i)=>[k,i]))};
 const poses=data.poses.map((row,i)=>{
  const pose={...row};
  for(const[field,namesKey]of [['scores','score_keys'],['percentiles','percentile_keys']])if(c.map_presence[field][i]){const view=mapView(c,field,c[namesKey],indexes[field],i);Object.defineProperty(pose,field,{enumerable:true,configurable:true,get:()=>view});}
  if(c.map_presence.contacts[i])Object.defineProperty(pose,'contacts',{enumerable:true,configurable:true,get:()=>c.contacts[i].map(j=>c.contact_vocab[j])});
  return pose;
 });
 const {columns,format,...envelope}=data;return{...envelope,poses};
}
return Object.freeze({validate,adapt});
});
