/* Low-level Mol* roots: no automatic presets and no inter-structure bonds. */
(function(root){'use strict';
const crypto={subtle:{digest:(algorithm,bytes)=>root.PoseHash.digest(algorithm,bytes)}};
const expr=(name,args)=>({head:{name},...(args===undefined?{}:{args})});
const prop=name=>expr('structure-query.atom-property.macromolecular.'+name);
const eq=(a,b)=>expr('core.rel.eq',[a,b]);
const or=items=>items.length===1?items[0]:expr('core.logic.or',items);
const and=items=>items.length===1?items[0]:expr('core.logic.and',items);
const residues=names=>expr('structure-query.generator.atom-groups',{'residue-test':or(names.map(n=>eq(prop('label_comp_id'),n)))});
class PoseViewer{
 constructor(host,status){this.host=host;this.element=typeof host==='string'?document.getElementById(host):host;this.status=status;this.instance=null;this.roots=new Map();this.version=0;this.reference=null;this.fetchController=null;this.pending=null;this.active=null;this.latest=null;this.initFlight=null;this.readiness='idle';this.initAttempts=0;this.lastError=null;this.hidden=false;this.geometryDirty=false;this.cameraReady=false;}
 init({retry=false}={}){
  if(this.readiness==='ready')return Promise.resolve(this);
  if(this.initFlight)return this.initFlight;
  if(this.readiness==='failed'&&!retry)return Promise.reject(this.lastError);
  if(this.initAttempts>=3)return Promise.reject(this.lastError||new Error('Viewer retry limit reached; reload the page.'));
  this.initAttempts++;this.readiness='warming';this.showStatus('Loading viewer…');
  this.initFlight=(async()=>{try{
   const vendor=await root.ViewerVendor.warm({retry});
   this.readiness='initializing';this.showStatus('Starting viewer…');
   this.instance=await vendor.Viewer.create(this.host,{layoutIsExpanded:false,layoutShowControls:false,layoutShowRemoteState:false,layoutShowSequence:false,layoutShowLog:false,layoutShowLeftPanel:false,collapseLeftPanel:true,viewportShowExpand:false,viewportShowSelectionMode:false,viewportShowAnimation:false,volumeStreamingDisabled:true});
   this.plugin.canvas3d?.setProps({renderer:{backgroundColor:0xffffff},marking:{enabled:false}});
   // Mol* keeps an unbounded diagnostic log even when its log panel is hidden.
   // Keep recent render diagnostics; scientific warnings remain in row data.
   this.logSubscription=this.plugin.events.log.subscribe(()=>{const log=this.plugin.log;if(log.entries.size>64)log.entries=log.entries.takeLast(64);});
   this.readiness='ready';this.lastError=null;this.showStatus('');this.drain();return this;
  }catch(error){this.readiness='failed';this.lastError=error;this.hide();this.failPending(error);this.showStatus('Viewer unavailable · '+error.message,true);throw error;}
  finally{this.initFlight=null;}})();return this.initFlight;
 }
 showStatus(text,retry=false){this.status.replaceChildren(document.createTextNode(text));if(retry&&this.initAttempts<3){const button=document.createElement('button');button.type='button';button.textContent='Retry viewer';button.style.pointerEvents='auto';button.onclick=()=>this.retry().catch(()=>{});this.status.append(' ',button);}}
 retry(){if(this.readiness==='ready'){if(!this.latest)return Promise.resolve(this);return this.set(this.latest.poses,this.latest.options);}if(this.initAttempts>=3)return Promise.reject(this.lastError||new Error('Viewer retry limit reached; reload the page.'));if(this.latest&&!this.pending)this.enqueue(this.latest.poses,this.latest.options).catch(()=>{});return this.init({retry:true});}
 failPending(error){if(this.pending){this.pending.reject(error);this.pending=null;}}
 hide(){this.hidden=true;this.element.style.visibility='hidden';}
 reveal(){this.hidden=false;this.element.style.visibility='visible';}
 get diagnostics(){const cells=this.plugin?[...this.plugin.state.data.cells.values()]:[];return Object.freeze({readiness:this.readiness,queuedCount:this.pending?1:0,activeCount:this.active?1:0,roots:Object.freeze([...this.roots.values()].map(item=>Object.freeze({id:item.pose.id,sha256:item.pose.structure_sha256,ref:item.ref}))),representationCount:cells.filter(c=>c.obj?.type?.typeClass==='Representation3D').length,receiptCount:root.cypCoordinateReceipts?.length||0,hidden:this.hidden});}
 get plugin(){return this.instance?.plugin;}
 identity(pose){return JSON.stringify([pose.structure_url,pose.structure_sha256,pose.source_state,pose.structure_format,pose.structure_key,pose.reference_role]);}
 compatible(poses){return poses.length===this.roots.size&&poses.every(p=>this.roots.get(p.id)?.identity===this.identity(p));}
 currentPose(pose){return this.latest?.poses.some(p=>p.id===pose.id&&this.identity(p)===this.identity(pose));}
 enqueue(poses,options){const version=++this.version;this.pending?.resolve({superseded:true});this.active?.resolve({superseded:true});return new Promise((resolve,reject)=>{this.pending={poses,options,version,resolve,reject};});}
 set(poses,options){
  // Copy the option snapshot: callers can change controls while Mol* is warming.
  poses=poses.slice();options={...options,colors:{...options.colors}};this.latest={poses,options};
  if(!this.compatible(poses))this.hide();
  if(this.fetchPose&&!this.currentPose(this.fetchPose))this.fetchController?.abort();
  const result=this.enqueue(poses,options);
  if(this.readiness==='failed'){this.failPending(this.lastError);this.showStatus('Viewer unavailable · '+this.lastError.message,true);}
  else this.drain();return result;
 }
 async drain(){if(this.active||!this.pending||this.readiness!=='ready')return;const job=this.pending;this.pending=null;this.active=job;
  try{await this.update(job.poses,job.options,job.version);job.resolve({superseded:job.version!==this.version});}
  catch(error){job.reject(error);}finally{this.active=null;this.drain();}
 }
 invalidate(){this.version++;this.latest=null;this.pending?.resolve({superseded:true});this.pending=null;this.active?.resolve({superseded:true});this.fetchController?.abort();this.hide();this.showStatus('Loading structures…');}
 cancel(){this.invalidate();this.showStatus('Loading stopped');}
 async remove(ref){await this.plugin.build().delete(ref).commit();}
 async component(structure,expression,key,type,color){const b=this.plugin.builders.structure;const component=await b.tryCreateComponentFromExpression(structure,expression,key,{label:key});const representation=component?await b.representation.addRepresentation(component,{type,color:'uniform',colorParams:{value:color},typeParams:type==='ball-and-stick'?{visuals:['element-sphere','intra-bond'],aromaticBonds:false,sizeFactor:.22}:type==='spacefill'?{sizeFactor:.38,alpha:1}:{alpha:.55}}):null;return{component,representation,color,visible:true};}
 async load(pose){const b=this.plugin.builders;const controller=new AbortController();this.fetchController=controller;this.fetchPose=pose;let text;try{const response=await fetch(PoseHelpers.localURL(pose.structure_url),{signal:controller.signal});if(!response.ok)throw new Error(`Structure ${response.status}: ${pose.id}`);let bytes=new Uint8Array(await response.arrayBuffer());if(bytes[0]===31&&bytes[1]===139){bytes=new Uint8Array(await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer());}const actual=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes))).map(v=>v.toString(16).padStart(2,'0')).join('');if(actual!==pose.structure_sha256)throw new Error('Decoded coordinate hash mismatch: '+pose.id);if(!this.currentPose(pose))return null;text=new TextDecoder('utf-8',{fatal:true}).decode(bytes);root.cypCoordinateReceipts??=[];root.cypCoordinateReceipts.push({id:pose.id,state:pose.source_state,expected:pose.structure_sha256,decoded_sha256:actual,bytes:bytes.length,url:pose.structure_url});root.cypCoordinateReceipts.splice(0,Math.max(0,root.cypCoordinateReceipts.length-16));}finally{if(this.fetchController===controller){this.fetchController=null;this.fetchPose=null;}}if(!this.currentPose(pose))return null;const raw=await b.data.rawData({data:text,label:pose.id});try{if(!this.currentPose(pose)){await this.remove(raw.ref);return null;}const format=pose.structure_format||(/\.cif(?:\?|$)/i.test(pose.structure_url)||text.trimStart().startsWith('data_')?'mmcif':'pdb');const trajectory=await b.structure.parseTrajectory(raw,format);if(!this.currentPose(pose)){await this.remove(raw.ref);return null;}const model=await b.structure.createModel(trajectory);const structure=await b.structure.createStructure(model);if(!this.currentPose(pose)){await this.remove(raw.ref);return null;}return{ref:raw.ref,structure,pose,identity:this.identity(pose),layers:{}};}catch(error){await this.remove(raw.ref);throw error;}}
 visibility(layer,on){if(layer.visible===on)return;for(const selector of [layer.component,layer.representation])if(selector)this.plugin.state.data.updateCellState(selector.ref,{isHidden:!on});layer.visible=on;}
 async layers(item,options){const {pose,structure,layers}=item;const color=parseInt((options.colors[pose.id]||'#547c8d').slice(1),16);const names={Ligand:pose.ligand_resnames||[pose.ligand_resname||'LIG'],Heme:pose.heme_resnames||[pose.heme_resname||'HEM'],Water:pose.water_resnames||[]};const enabled={Ligand:options.ligand,Heme:options.heme,Water:!!pose.water_resnames?.length,Reference:options.protein&&(pose.id===this.reference||!!pose.reference_role),Pocket:options.pocket};
 for(const [key,on]of Object.entries(enabled)){
  if(on&&!layers[key]){
   if(key==='Reference'){const component=await this.plugin.builders.structure.tryCreateComponentStatic(structure,'polymer','Reference');const representation=component?await this.plugin.builders.structure.representation.addRepresentation(component,{type:'cartoon',color:'uniform',colorParams:{value:pose.reference_role==='revealed'?0x83bdbe:0xb9bfc2},typeParams:{alpha:.5}}):null;layers[key]={component,representation,visible:true};}
   else if(key==='Pocket'){const tests=(pose.pocket_residues||[]).filter(r=>r.chain!==undefined&&Number.isInteger(r.residue)).map(r=>and([eq(prop('auth_asym_id'),r.chain),eq(prop('auth_seq_id'),r.residue)]));if(tests.length)layers[key]=await this.component(structure,expr('structure-query.generator.atom-groups',{'residue-test':or(tests)}),key,'ball-and-stick',0x9aabb0);}
   else layers[key]=await this.component(structure,residues(names[key]),key,key==='Water'?'spacefill':'ball-and-stick',key==='Ligand'?color:key==='Water'?0xe53935:0x926c56);
  }
  if(layers[key])this.visibility(layers[key],!!on);
 }
 const ligand=layers.Ligand;if(ligand?.representation&&ligand.color!==color){await this.plugin.build().to(ligand.representation).update(old=>({...old,colorTheme:{name:'uniform',params:{value:color}}})).commit();ligand.color=color;}
 }
 async settleCamera(canvas,version){const previous=canvas.notifyDidDraw;let drawn=false;canvas.notifyDidDraw=true;const subscription=canvas.didDraw.subscribe(()=>{drawn=true;});const started=performance.now();try{do{if(version!==this.version)return;canvas.requestDraw();await new Promise(requestAnimationFrame);if(performance.now()-started>10000)throw new Error('Viewer camera did not settle');}while(!drawn||canvas.camera.transition.inTransition||(canvas.commitQueueSize?.value||0)>0);}finally{subscription.unsubscribe();canvas.notifyDidDraw=previous;}}
 async update(poses,options,version){if(version!==this.version||!this.plugin)return;const p=this.plugin;p.canvas3d?.setProps({camera:{manualReset:!!options.lock},trackball:{autoAdjustMinMaxDistance:{name:'off',params:{}}}});const camera=options.lock?(this.cameraReady&&this.roots.size?p.canvas3d?.camera.getSnapshot():this.lockedCamera):null;if(camera)this.lockedCamera=camera;if(!options.lock)this.lockedCamera=null;if(!poses.some(x=>x.id===this.reference))this.reference=poses[0]?.id||null;const desired=new Map(poses.map(p=>[p.id,this.identity(p)]));const geometry=[...this.roots].some(([id,item])=>desired.get(id)!==item.identity)||poses.some(p=>!this.roots.has(p.id));const layerSignature=JSON.stringify([options.ligand,options.heme,options.protein,options.pocket,this.reference]);const structural=geometry||this.geometryDirty||this.layerSignature!==layerSignature;this.status.textContent=poses.length?'Loading structures…':'';let guard=0;
 const restore=()=>{if(camera)p.canvas3d?.camera.setState(camera,0);};if(camera&&structural)guard=setInterval(restore,35);
 try{const removed=[...this.roots].filter(([id,item])=>desired.get(id)!==item.identity);if(removed.length){const update=p.build();for(const[,item]of removed)update.delete(item.ref);await update.commit();for(const[id]of removed)this.roots.delete(id);this.geometryDirty=true;}
 let completed=0;for(const pose of poses){if(version!==this.version)return;completed++;const signature=JSON.stringify([options.ligand,options.heme,options.protein&&(pose.id===this.reference||!!pose.reference_role),options.pocket,options.colors[pose.id]]);if(this.roots.get(pose.id)?.layerOptions===signature)continue;this.status.textContent=`Loading ${completed} / ${poses.length} poses…`;await new Promise(resolve=>setTimeout(resolve,0));if(version!==this.version)return;if(!this.roots.has(pose.id)){const item=await this.load(pose);if(!item)return;if(!this.currentPose(pose)){await this.remove(item.ref);return;}this.roots.set(pose.id,item);this.geometryDirty=true;}if(version!==this.version)return;const item=this.roots.get(pose.id);await this.layers(item,options);item.layerOptions=signature;document.getElementById('viewerCount').textContent=`${this.roots.size} / ${poses.length} poses`;restore();}
 if(version!==this.version)return;if(structural)p.canvas3d?.requestResize();if(!camera&&poses.length&&structural){p.managers.camera.reset();await this.settleCamera(p.canvas3d,version);}if(version!==this.version)return;restore();this.layerSignature=layerSignature;this.geometryDirty=false;if(poses.length)this.cameraReady=true;this.reveal();this.status.textContent=options.pocket&&!poses.some(x=>x.pocket_residues?.length)?'Pocket residues unavailable':'';
 }catch(error){if(version===this.version){this.hide();this.status.textContent=error.message;throw error;}}finally{if(guard){await new Promise(r=>setTimeout(r,180));clearInterval(guard);restore();}if(version===this.version)document.getElementById('viewerCount').textContent=`${this.roots.size} poses`;}}
 fit(){this.plugin?.managers.camera.reset();}
 async screenshot(){const helper=this.plugin?.helpers.viewportScreenshot;if(!helper)throw new Error('Viewer image export unavailable');return helper.getImageDataUri();}
}
root.PoseViewer=PoseViewer;
})(globalThis);
