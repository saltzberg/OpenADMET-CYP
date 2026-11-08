/* Pure UI helpers. No scientific values or projections are generated here. */
(function(root,factory){const api=factory();if(typeof module==='object')module.exports=api;else root.PoseHelpers=api;})(typeof globalThis!=='undefined'?globalThis:this,function(){
'use strict';
const finite=x=>typeof x==='number'&&Number.isFinite(x);
const xy=x=>Array.isArray(x)&&x.length===2&&x.every(finite);
const fmt=x=>finite(x)?Number(x).toLocaleString('en-US',{maximumSignificantDigits:4}):'—';
const escape=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function element(doc,tag,text,attrs={}){const el=doc.createElement(tag);if(text!==undefined)el.textContent=String(text);for(const[k,v]of Object.entries(attrs))el.setAttribute(k,String(v));return el;}
function apiURL(path,params={}){const q=new URLSearchParams();for(const[k,v]of Object.entries(params))if(v!==undefined&&v!==null)q.set(k,String(v));return path+(q.size?'?'+q:'');}
function localURL(value){if(typeof value!=='string'||!value.startsWith('/')||value.startsWith('//')||value.includes('\\'))throw new Error('Expected a local asset URL');return value;}
function csv(rows){const keys=[...new Set(rows.flatMap(r=>Object.keys(r)))];const cell=v=>'"'+String(v==null?'':typeof v==='object'?JSON.stringify(v):v).replace(/"/g,'""')+'"';return[keys.map(cell).join(','),...rows.map(r=>keys.map(k=>cell(r[k])).join(','))].join('\r\n');}
function clusterDisplay(cluster,poses){const representative=poses.find(p=>p.id===cluster?.representative_id);return representative&&representative.available!==false&&representative.structure_url?[representative.id]:[];}
function sorted(rows,key,direction=1){return [...rows].sort((a,b)=>{const x=key(a),y=key(b);if(x==null)return y==null?0:1;if(y==null)return-1;return direction*(typeof x==='number'&&typeof y==='number'?x-y:String(x).localeCompare(String(y),undefined,{numeric:true}));});}
function copyPose(p,updates={}){const descriptors=Object.getOwnPropertyDescriptors(p);for(const[key,value]of Object.entries(updates))descriptors[key]={value,enumerable:true,writable:true,configurable:true};return Object.defineProperties({},descriptors);}
function checkReference(ref,data){
 if(ref==null)return;
 if(ref.role!=='revealed_reference'||ref.snapshot!==data.snapshot||ref.source_snapshot!==data.source_snapshot||ref.ligand!==data.ligand||ref.state!==data.state||!/^[0-9a-f]{64}$/.test(ref.revision||'')||!/^[0-9a-f]{64}$/.test(ref.source_pdb_sha256||''))throw new Error('Mixed or invalid revealed reference');
 for(const field of ['pose','receptor','energy']){const p=ref[field];if(!p||!['AVAILABLE','UNAVAILABLE'].includes(p.status)||(p.status==='AVAILABLE'?!xy(p.xy):p.xy!==null))throw new Error('Invalid revealed reference coordinates');}
 if(ref.pose.coordinate_state!==data.state)throw new Error('Mixed reference coordinate state');
}
function referencePoint(data,request,plot){
 const ref=data?.revealed_reference;
 if(!ref||ref.ligand!==request.ligand||ref.state!==request.state)return null;
 const field={posePlot:'pose',receptorPlot:'receptor',energyPlot:'energy'}[plot],point=ref[field];
 if(field==='pose'&&(point.metric!==request.metric||(request.metric==='selected_atoms'&&JSON.stringify(point.atom_indices)!==JSON.stringify(request.atoms))))return{...point,status:'PENDING',xy:null,reason:'Reference projection pending for current metric / atoms'};
 return point;
}
function applyAnalysis(data,result){checkReference(result.revealed_reference,{snapshot:result.snapshot,source_snapshot:result.source_snapshot,ligand:data.ligand,state:data.state});if(!Array.isArray(result.pose_ids)||!Array.isArray(result.xy)||!Array.isArray(result.labels)||result.xy.length!==result.pose_ids.length||result.labels.length!==result.pose_ids.length)throw new Error('Invalid analysis shape');let spec=result.analysis_spec;if(spec?.pose_ids_source!==undefined){if(spec.pose_ids_source!=='pose_ids')throw new Error('Unknown analysis identity reference');const {pose_ids_source,...fields}=spec;spec={...fields,pose_ids:result.pose_ids};}const updates=new Map(result.pose_ids.map((id,i)=>[id,{pose_xy:xy(result.xy[i])?result.xy[i]:null,cluster:result.labels[i]}]));return{...data,poses:data.poses.map(p=>copyPose(p,{pose_xy:null,cluster:null,...updates.get(p.id)})),clusters:result.clusters||[],projection:result.projection,warnings:result.warnings||[],analysis_spec:spec,revealed_reference:result.revealed_reference??null};}
function checkLigand(data){if(!data||!Array.isArray(data.poses)||!Array.isArray(data.atoms)||!data.coverage)throw new Error('Invalid ligand response');checkReference(data.revealed_reference,data);const ids=new Set();for(const p of data.poses){if(typeof p.id!=='string'||ids.has(p.id))throw new Error('Missing or duplicate pose ID');ids.add(p.id);}return data;}
class Latest{constructor(){this.version=0;this.controller=null;}cancel(){this.controller?.abort();this.version++;this.controller=null;}start(){this.cancel();this.controller=new AbortController();return{version:this.version,signal:this.controller.signal};}current(ticket){return ticket.version===this.version&&!ticket.signal.aborted;}}
function atomIndices(indices,atoms){const allowed=new Set(atoms.map(a=>a.index));return [...new Set(indices.filter(i=>Number.isInteger(i)&&allowed.has(i)))].sort((a,b)=>a-b);}
function pinnedURL(value,session){
 if(!session?.api_base||typeof value!=='string')throw new Error('Snapshot not initialized');
 const prefix='/articles/admet-cyp-structure/cofolding/';let url=value;
 if(url.startsWith('./api/'))url='/'+url.slice(2);else if(url.startsWith('api/'))url='/'+url;
 localURL(url);if(decodeURIComponent(url.split('?')[0]).split('/').includes('..'))throw new Error('Invalid asset path');
 if(url.startsWith(session.api_base+'/'))return url;
 if(url.startsWith(prefix+'cyp-main/'))return url;
 if(url.startsWith(prefix+'api/'))url='/'+url.slice(prefix.length);
 if(url.startsWith('/api/v1/')){
  const short='/api/v1/'+session.snapshot+'/';
  if(!url.startsWith(short))throw new Error('Snapshot mismatch — Reload required');
  return session.api_base+'/'+url.slice(short.length);
 }
 if(!url.startsWith('/api/'))throw new Error('Expected API asset');return session.api_base+url.slice(4);
}
function adaptLigand(data,session){return{...data,poses:data.poses.map(p=>{
 const pose=copyPose(p,p.structure_url?{structure_url:pinnedURL(p.structure_url,session)}:{});
 if(!Object.hasOwn(pose,'_searchText'))Object.defineProperty(pose,'_searchText',{value:[pose.id,pose.method,...(pose.contacts||[])].join(' ').toLowerCase(),enumerable:false});
 return pose;
})};}
function analysisMatches(data,request,session){
 const a=data?.analysis_spec;if(!a||a.snapshot!==session.snapshot||data.snapshot!==session.snapshot||a.ligand!==request.ligand||a.state!==request.state||a.metric!==request.metric||a.threshold!==request.threshold)return false;
 if(a.source_snapshot!==(session.analysis_source_snapshot??session.source_snapshot)||a.kernel_sha256!==session.analysis_kernel_sha256||a.projection_version!==session.projection_version||a.cluster_algorithm!==session.cluster_algorithm)return false;
 if(!/^[0-9a-f]{64}$/.test(a.population_sha256||''))return false;
 const atoms=request.metric==='selected_atoms'?request.atoms:[];
 if(!Array.isArray(a.atom_indices)||a.atom_indices.length!==atoms.length||a.atom_indices.some((x,i)=>x!==atoms[i]))return false;
 return Array.isArray(a.pose_ids)&&a.pose_ids.length===data.poses.length&&a.pose_ids.every((id,i)=>id===data.poses[i].id);
}
function debounce(callback,delay=150){let timer=null,args=[];const invoke=()=>{timer=null;callback(...args);};const run=(...next)=>{args=next;clearTimeout(timer);timer=setTimeout(invoke,delay);};run.flush=()=>{if(timer!==null){clearTimeout(timer);invoke();}};run.cancel=()=>{clearTimeout(timer);timer=null;};return run;}
return{finite,xy,fmt,escape,element,apiURL,localURL,csv,clusterDisplay,sorted,applyAnalysis,checkLigand,checkReference,referencePoint,Latest,atomIndices,copyPose,pinnedURL,adaptLigand,analysisMatches,debounce};
});
