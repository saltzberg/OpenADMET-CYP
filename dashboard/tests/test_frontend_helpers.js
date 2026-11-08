/* Functional helper tests. Synthetic scalar controls plus existing pose metadata. */
'use strict';
const test=require('node:test'),assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
const H=require('../static/helpers.js');
const session={snapshot:'build1',source_snapshot:'science1',api_base:'/articles/admet-cyp-structure/cofolding/api/v1/build1',analysis_kernel_sha256:'a'.repeat(64),projection_version:'pcoa-canonical-eigenspaces-v2',cluster_algorithm:'average-linkage-original-distances-v1'};
function freeze(x){if(x&&typeof x==='object'){Object.values(x).forEach(freeze);Object.freeze(x);}return x;}
test('pinning preserves semantic queries and is idempotent for supported URL forms',()=>{
 for(const url of ['/api/ligand?ligand=OCNT-2312554&state=paired',session.api_base+'/ligand?state=raw','api/ligand?state=relaxed','./api/ligand','/articles/admet-cyp-structure/cofolding/api/ligand','/api/v1/build1/ligand','/articles/admet-cyp-structure/cofolding/cyp-main/raw/x.cif.gz']){
  const once=H.pinnedURL(url,session);assert.equal(H.pinnedURL(once,session),once);
 }
 assert.equal(H.pinnedURL('/api/ligand?x=a%26b&state=raw',session),session.api_base+'/ligand?x=a%26b&state=raw');
 for(const url of ['https://external.invalid/api','//external.invalid/api','/api/v1/old/ligand','/articles/admet-cyp-structure/cofolding/api/v1/old/ligand','/api/../private','/api/%2e%2e/private'])assert.throws(()=>H.pinnedURL(url,session));
});
test('existing pose adaptation and analysis do not mutate any nested inputs',()=>{
 const raw=JSON.parse(fs.readFileSync(path.join(process.env.CYP_TEST_APP || path.resolve(__dirname,'..'),'data/features/OCNT-2312554/raw.json')))[0];
 const pose={...raw,structure_url:'/api/structure/example',scores:{control:.12345678901234567},percentiles:{control:null},contacts:raw.contacts||[],pose_xy:[1,2],cluster:1};
 const data=freeze({poses:[pose],atoms:[],coverage:{shown:1},snapshot:session.snapshot});const before=JSON.stringify(data);
 const adapted=H.adaptLigand(data,session);assert.notEqual(adapted.poses[0],pose);assert.equal(JSON.stringify(data),before);assert.equal(adapted.poses[0]._searchText,[pose.id,pose.method,...pose.contacts].join(' ').toLowerCase());assert(!JSON.stringify(adapted).includes('_searchText'));
 const analysis=freeze({pose_ids:[pose.id],xy:[[3,4]],labels:[2],clusters:[],projection:{version:'test'}});
 const result=H.applyAnalysis(freeze(adapted),analysis);assert.equal(JSON.stringify(data),before);assert.deepEqual(result.poses[0].pose_xy,[3,4]);assert.equal(result.poses[0].scores.control,pose.scores.control);assert.equal(result.poses[0]._searchText,adapted.poses[0]._searchText);
 const reference=freeze({...analysis,analysis_spec:{snapshot:session.snapshot,pose_ids_source:'pose_ids'}});const refBefore=JSON.stringify(reference);const resolved=H.applyAnalysis(adapted,reference);assert.deepEqual(resolved.analysis_spec,{snapshot:session.snapshot,pose_ids:[pose.id]});assert.equal(JSON.stringify(reference),refBefore);assert.throws(()=>H.applyAnalysis(adapted,{...analysis,analysis_spec:{pose_ids_source:'unknown'}}));
});
test('full analysis specification must match every requested scientific identity',()=>{
 const request={ligand:'OCNT-2312554',state:'relaxed',metric:'rmsd',atoms:[],threshold:2};
 const spec={snapshot:session.snapshot,source_snapshot:session.source_snapshot,ligand:request.ligand,state:request.state,metric:request.metric,atom_indices:[],pose_ids:['a','b'],threshold:2,kernel_sha256:session.analysis_kernel_sha256,projection_version:session.projection_version,cluster_algorithm:session.cluster_algorithm,population_sha256:'b'.repeat(64)};
 const data={snapshot:session.snapshot,poses:[{id:'a'},{id:'b'}],analysis_spec:spec};assert(H.analysisMatches(data,request,session));
 assert(H.analysisMatches(data,request,{...session,analysis_source_snapshot:session.source_snapshot,scientific_source:{snapshot:'different-logical-display-id'}}));
 assert(!H.analysisMatches(data,request,{...session,analysis_source_snapshot:'different-logical-display-id'}));
 for(const field of ['snapshot','source_snapshot','ligand','state','metric','kernel_sha256','projection_version','cluster_algorithm','population_sha256'])assert(!H.analysisMatches({...data,analysis_spec:{...spec,[field]:'mismatch'}},request,session),field);
 for(const update of [{threshold:2.25},{atom_indices:[0]},{pose_ids:['b','a']},{pose_ids:['a']},{pose_ids:['a','c']}])assert(!H.analysisMatches({...data,analysis_spec:{...spec,...update}},request,session));
 assert(!H.analysisMatches(data,{...request,metric:'contacts'},session));
 assert(!H.analysisMatches(data,{...request,threshold:2.1},session));
});
test('150 ms input burst coalesces; export flush is synchronous; cancel removes pending work',async()=>{
 let calls=[];const run=H.debounce(x=>calls.push(x),150);run('a');run('ab');run('abc');assert.deepEqual(calls,[]);await new Promise(r=>setTimeout(r,180));assert.deepEqual(calls,['abc']);run('export');run.flush();assert.deepEqual(calls,['abc','export']);await new Promise(r=>setTimeout(r,180));assert.equal(calls.length,2);run('cancel');run.cancel();await new Promise(r=>setTimeout(r,180));assert.equal(calls.length,2);
});
