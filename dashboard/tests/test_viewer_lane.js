/* Unit contracts only. Real coordinates/WebGL are checked by the Python browser suite. */
'use strict';
const test=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const vm=require('node:vm');
const path=require('node:path');
const staticPath=path.resolve(__dirname,'../static');
function context(){
 const timers=new Map(),scripts=[];
 function element(){return {style:{},children:[],replaceChildren(...nodes){this.children=nodes;},append(...nodes){this.children.push(...nodes);},remove(){this.removed=true;}};}
 const host=element(),status=element();
 const scope={console,AbortController,setTimeout(fn){const id={};timers.set(id,fn);return id;},clearTimeout(id){timers.delete(id);},document:{getElementById(id){return id==='viewer'?host:status;},createElement:element,createTextNode:text=>({text}),head:{append(script){scripts.push(script);}}}};
 scope.globalThis=scope;vm.createContext(scope);
 vm.runInContext(fs.readFileSync(path.join(staticPath,'vendor-loader.js'),'utf8'),scope);
 return {scope,timers,scripts,host,status};
}
const tick=()=>new Promise(resolve=>setImmediate(resolve));
test('UMD loader is single-flight, times out, and stops after three attempts',async()=>{
 const {scope:s,timers,scripts}=context();
 const first=s.ViewerVendor.warm();assert.equal(first,s.ViewerVendor.warm());assert.equal(scripts.length,1);
 let failed=assert.rejects(first,/timed out/);[...timers.values()][0]();await failed;
 assert.equal(s.ViewerVendor.state,'failed');assert.equal(scripts[0].removed,true);
 assert.equal(s.ViewerVendor.warm(),first);
 for(let attempt=2;attempt<=3;attempt++){
  const next=s.ViewerVendor.warm({retry:true});failed=assert.rejects(next,/could not load/);scripts.at(-1).onerror();await failed;
 }
 await assert.rejects(s.ViewerVendor.warm({retry:true}));assert.equal(scripts.length,3);
 assert.equal(s.ViewerVendor.attempts,3);
});
test('late successful UMD execution after timeout can be used by an explicit retry',async()=>{
 const {scope:s,timers}=context();const first=s.ViewerVendor.warm();const failed=assert.rejects(first);[...timers.values()][0]();await failed;
 s.molstar={Viewer:{}};assert.equal(await s.ViewerVendor.warm({retry:true}),s.molstar);
});
test('Plotly readiness rejects a missing global and returns the existing global',async()=>{
 const {scope:s}=context();await assert.rejects(s.PlotlyVendor.ready(),/unavailable/);s.Plotly={};assert.equal(await s.PlotlyVendor.ready(),s.Plotly);
});
test('viewer retains latest pre-init request; each superseded caller settles',async()=>{
 const {scope:s,host,status}=context();vm.runInContext(fs.readFileSync(path.join(staticPath,'viewer.js'),'utf8'),s);
 let release,logCallback;let creates=0;const plugin={canvas3d:{setProps(){}},state:{data:{cells:new Map()}},log:{entries:{size:100,takeLast(size){return{size}}}},events:{log:{subscribe(fn){logCallback=fn;return{unsubscribe(){}}}}}};
 s.ViewerVendor={warm:()=>Promise.resolve({Viewer:{create:async()=>{creates++;await new Promise(r=>release=r);return {plugin}}}})};
 const v=new s.PoseViewer(host,status);const updates=[];
 v.update=async(poses,options,version)=>{updates.push({poses,options,version})};
 const options={colors:{},lock:true};const a=v.set([{id:'unit-a'}],options);const b=v.set([{id:'unit-b'}],options);
 assert.equal((await a).superseded,true);assert.equal(v.diagnostics.queuedCount,1);
 const init=v.init();assert.equal(init,v.init());await tick();release();await init;await b;await tick();
 assert.equal(creates,1);assert.equal(updates.length,1);assert.equal(updates[0].poses[0].id,'unit-b');
 assert.equal(v.plugin,plugin);assert.equal(v.readiness,'ready');assert.equal(v.diagnostics.queuedCount,0);
 logCallback();assert.equal(plugin.log.entries.size,64);
});
test('viewer coalesces active work plus one pending request; cancel settles pending callers',async()=>{
 const {scope:s,host,status}=context();vm.runInContext(fs.readFileSync(path.join(staticPath,'viewer.js'),'utf8'),s);
 const v=new s.PoseViewer(host,status);v.readiness='ready';let release;let runs=0;
 v.update=async()=>{runs++;await new Promise(r=>release=r)};
 const a=v.set([{id:'unit-a'}],{colors:{}});const jobs=[];
 for(let i=0;i<100;i++)jobs.push(v.set([{id:'unit-'+i}],{colors:{}}));
 assert.equal(v.diagnostics.queuedCount,1);assert.equal(v.diagnostics.activeCount,1);assert.equal((await a).superseded,true);
 v.cancel();await Promise.all(jobs);release();await tick();assert.equal(runs,1);assert.equal(v.diagnostics.queuedCount,0);assert.equal(v.diagnostics.activeCount,0);assert.equal(host.style.visibility,'hidden');
});
test('viewer init failures retain latest request and bound explicit retry',async()=>{
 const {scope:s,host,status}=context();vm.runInContext(fs.readFileSync(path.join(staticPath,'viewer.js'),'utf8'),s);
 s.ViewerVendor={warm:async()=>{throw new Error('unit bundle failure')}};
 const v=new s.PoseViewer(host,status);const pending=v.set([{id:'unit-a'}],{colors:{}});const failed=assert.rejects(pending);
 await assert.rejects(v.init());await failed;
 for(let i=0;i<2;i++)await assert.rejects(v.retry());
 await assert.rejects(v.retry());assert.equal(v.initAttempts,3);assert.equal(v.diagnostics.queuedCount,0);assert.equal(v.latest.poses[0].id,'unit-a');
});
