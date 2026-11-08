/* Direct, chemically typed, angular hydrogen bonds; source-state-bound results. */
(function(root,factory){'use strict';const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;else root.ResidueContacts=api;})(typeof globalThis!=='undefined'?globalThis:this,function(){'use strict';
 const crypto={subtle:{digest:(algorithm,bytes)=>(globalThis.PoseHash||require('./sha256.js')).digest(algorithm,bytes)}};
 const definition='Direct sidechain hydrogen bond: chemically compatible donor/acceptor, donor–acceptor distance ≤3.5 Å and donor–H–acceptor angle ≥130°. SER119 OG; ARG105 NE/NH1/NH2 as donors. No backbone or water bridge. Uses the exact row coordinate state and recorded hydrogen preparation. Missing evaluation is Unknown.';
 const metrics=Object.freeze(['ser119_contact','arg105_contact'].map((id,i)=>Object.freeze({id,label:['SER119','ARG105'][i]+' direct H-bond',kind:'contact',encoding:'raw',range:[0,1],definition,group:'Direct hydrogen bonds'})));
 const base='/articles/admet-cyp-structure/cofolding/static/residue-hbonds/';
 const states=['raw','relaxed'],listeners=new Set(),cache=new Map();
 let manifest=null,active=null,controller=null,generation=0,pinnedSession=null,pinnedIdentity=null,pinnedSource=null;
 let status=Object.freeze({state:'idle',loading:false,ready:false,failed:false,ligand:null,revision:null,error:null});
 const freeze=data=>{if(data&&typeof data==='object'&&!Object.isFrozen(data)){Object.values(data).forEach(freeze);Object.freeze(data);}return data;};
 function emit(state,ligand,revision,error=null){status=Object.freeze({state,loading:state==='loading',ready:state==='ready',failed:state==='failed',ligand,revision,error});for(const fn of listeners){try{fn(status);}catch(e){console.error('H-bond subscriber failed',e);}}}
 function subscribe(fn){if(typeof fn!=='function')throw new TypeError('Expected H-bond callback');listeners.add(fn);return()=>listeners.delete(fn);}
 function invalidate(){generation++;controller?.abort();controller=null;active=null;emit('idle',null,null);}
 const snapshot=value=>typeof value==='string'?value:value?.snapshot;
 function identity(session){return JSON.stringify([session?.snapshot,snapshot(session?.source_snapshot??session?.bootstrap?.source_snapshot),session?.source_scientific_revision??session?.bootstrap?.source_scientific_revision]);}
 function compatible(session,m){
  const source=snapshot(session?.source_snapshot??session?.bootstrap?.source_snapshot);
  const revision=session?.source_scientific_revision??session?.bootstrap?.source_scientific_revision;
  if(source)return source===m.source_snapshot&&(!revision||revision===m.source_scientific_revision);
  if(revision)return revision===m.source_scientific_revision;
  return session?.snapshot===m.original_application_snapshot||session?.snapshot===m.source_snapshot||session?.snapshot===m.source_scientific_revision;
 }
 const sameStates=value=>Array.isArray(value)&&value.length===2&&states.every((state,i)=>value[i]===state);
 const sourceIdentity=m=>JSON.stringify([m.source_index_sha256,m.source_snapshot,m.source_scientific_revision,m.original_application_snapshot]);
 function validateManifest(m){
  if(m?.schema!=='residue-hbonds-manifest'||m.schema_version!==1||!sameStates(m.states)||!m.shards||!m.source_metadata||!m.source_index_sha256||!m.source_snapshot||!m.source_scientific_revision||Object.keys(m.shards).length!==m.dimensions?.ligands)throw new Error('Invalid H-bond manifest schema');
  for(const [ligand,entry] of Object.entries(m.shards))if(!/^[a-f0-9]{64}$/.test(entry.sha256)||entry.file!==ligand+'.'+entry.sha256+'.json'||!Number.isInteger(entry.bytes)||entry.bytes<=0||entry.dimensions?.states!==2)throw new Error('Invalid H-bond shard descriptor');
  return freeze(m);
 }
 async function getText(url,signal){const response=await fetch(url,{cache:'no-cache',signal});if(!response.ok)throw new Error('Hydrogen-bond data unavailable: HTTP '+response.status);return response.text();}
 async function sha256(bytes){return Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes)),x=>x.toString(16).padStart(2,'0')).join('');}
 function validateShard(shard,m,ligand,entry){
  if(shard?.schema!=='residue-hbonds-shard'||shard.schema_version!==1||shard.ligand!==ligand||!sameStates(shard.states)||!shard.poses||Object.keys(shard.poses).length!==entry.dimensions.rows||shard.dimensions?.rows!==entry.dimensions.rows||shard.dimensions?.states!==2)throw new Error('Invalid H-bond shard schema/dimensions');
  for(const field of ['source_index_sha256','source_snapshot','source_scientific_revision','original_application_snapshot'])if(shard[field]!==m[field])throw new Error('H-bond shard source mismatch — Reload required');
  const counts={raw:0,relaxed:0};
  for(const [key,row] of Object.entries(shard.poses)){
   if(!states.includes(row.source_state)||!key.startsWith(row.source_state+'|')||!(/^[a-f0-9]{64}$/).test(row.structure_sha256)||!metrics.every(({id})=>row[id]===null||row[id]===0||row[id]===1))throw new Error('Invalid H-bond row');
   counts[row.source_state]++;
  }
  for(const state of states)if(counts[state]!==entry.state_coverage[state].rows||counts[state]!==shard.state_coverage?.[state]?.rows)throw new Error('H-bond state coverage mismatch');
  return freeze(shard);
 }
 async function load(session,ligand){
  const ticket=++generation;controller?.abort();controller=new AbortController();active=null;
  const signal=controller.signal,binding=identity(session),current=()=>ticket===generation&&binding===identity(session);
  emit('loading',ligand,session?.snapshot??null);
  try{
   if(!session?.snapshot||typeof ligand!=='string')throw new Error('H-bond load requires a session and ligand');
   // The mutable manifest is revalidated on every explicit load, including retries and LRU hits.
   const m=validateManifest(JSON.parse(await getText(base+'manifest.json',signal)));
   if(!current())throw new DOMException('Superseded H-bond request','AbortError');
   if(!compatible(session,m)){cache.clear();manifest=null;throw new Error('H-bond source snapshot mismatch — Reload required');}
   if(pinnedSession===session&&pinnedIdentity===binding&&pinnedSource!==sourceIdentity(m)){cache.clear();manifest=null;throw new Error('H-bond source revision changed in this session — Reload required');}
   pinnedSession=session;pinnedIdentity=binding;pinnedSource=sourceIdentity(m);
   if(manifest?.source_index_sha256!==m.source_index_sha256)cache.clear();
   manifest=m;
   const entry=m.shards[ligand];if(!entry)throw new Error('No H-bond shard for ligand '+ligand);
   let cached=cache.get(ligand);
   if(cached?.sha256!==entry.sha256){
    // Evict before parsing a new shard: retained row data never exceeds two ligands.
    cache.delete(ligand);while(cache.size>=2)cache.delete(cache.keys().next().value);
    const text=await getText(base+entry.file,signal),bytes=new TextEncoder().encode(text);
    if(bytes.byteLength!==entry.bytes||await sha256(bytes)!==entry.sha256)throw new Error('H-bond shard SHA-256 mismatch');
    if(!current())throw new DOMException('Superseded H-bond request','AbortError');
    cached={sha256:entry.sha256,bytes:entry.bytes,shard:validateShard(JSON.parse(text),m,ligand,entry)};
   }
   if(!current())throw new DOMException('Superseded H-bond request','AbortError');
   cache.delete(ligand);cache.set(ligand,cached);
   active={session,binding,ligand,shard:cached.shard};
   emit('ready',ligand,session.snapshot);
   return Object.freeze({manifest:m,shard:cached.shard});
  }catch(e){if(ticket===generation){active=null;if(e.name!=='AbortError')emit('failed',ligand,session?.snapshot??null,e.message);else emit('idle',null,null);}throw e;}
 }
 function value(pose,id){
  if(!metrics.some(m=>m.id===id)||!pose||pose.available===false||!pose.structure_sha256||!['raw','relaxed'].includes(pose.source_state))return null;
  const state=pose.source_state,suffix='::'+state,pid=pose.id?.endsWith(suffix)?pose.id.slice(0,-suffix.length):pose.id;
  if(!active||status.state!=='ready'||identity(active.session)!==active.binding)return null;
  const row=active.shard.poses[state+'|'+pid];
  if(!row||row.structure_sha256!==pose.structure_sha256)return null;
  const result=row[id];return result===0||result===1?result:null;
 }
 const label=value=>value===1?'H-bond':value===0?'No H-bond':'Unknown';
 const color=value=>value===1?'#0072b2':value===0?'#d89000':'#bdc2c5';
 const matches=(pose,id,filter)=>filter===''||filter==null||(filter==='unknown'?value(pose,id)===null:filter==='1'?value(pose,id)===1:filter==='0'?value(pose,id)===0:false);
 return Object.freeze({metrics,definition,value,label,color,matches,load,invalidate,subscribe,get status(){return status;},get cacheInfo(){return Object.freeze({ligands:Object.freeze([...cache.keys()]),shards:cache.size,totalBytes:[...cache.values()].reduce((n,v)=>n+v.bytes,0)});},get provenance(){return manifest?Object.freeze({...manifest.source_metadata,source_index_sha256:manifest.source_index_sha256,source_snapshot:manifest.source_snapshot,source_scientific_revision:manifest.source_scientific_revision,state_coverage:manifest.state_coverage,status}):null;}});
});
