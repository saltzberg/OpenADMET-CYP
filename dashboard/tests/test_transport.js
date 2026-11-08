'use strict';
const test=require('node:test'),assert=require('node:assert/strict'),crypto=require('node:crypto'),{spawnSync}=require('node:child_process'),path=require('node:path');
const T=require('../static/transport.js'),H=require('../static/helpers.js');
test('served UI bundle matches both canonical modules',()=>{const run=spawnSync((process.env.PYTHON || 'python3'),[path.resolve(__dirname,'../build_assets.py'),'--check'],{encoding:'utf8'});assert.equal(run.status,0,run.stdout+run.stderr);});
function wire(){const ids=['real-id-a','real-id-b'];return{format:'columns-v1',snapshot:'s',ligand:'l',state:'paired',poses:ids.map(id=>({id})),columns:{version:1,snapshot:'s',ligand:'l',state:'paired',pose_ids:ids,population_sha256:crypto.createHash('sha256').update(JSON.stringify(ids)).digest('hex'),dimensions:{poses:2,scores:2,percentiles:1,contacts:2},score_keys:['value','optional'],percentile_keys:['value'],scores:[[.12345678901234567,null],[null,9]],percentiles:[[null,.5]],presence:{scores:[null,[false,true]],percentiles:[null]},map_presence:{scores:[true,true],percentiles:[true,false],contacts:[true,true]},contact_vocab:['second','first'],contacts:[[1,0,1],[]]}};}
test('synchronous map access preserves null, absent keys, precision and ordered duplicate contacts',async()=>{
 const source=wire(),before=JSON.stringify(source),d=await T.adapt(source);assert.equal(d.poses[0].scores.value,source.columns.scores[0][0]);assert.equal(d.poses[1].scores.value,null);assert(!('optional' in d.poses[0].scores));assert(!('percentiles' in d.poses[1]));assert.deepEqual(d.poses[0].contacts,['first','second','first']);assert.deepEqual(Object.entries(d.poses[0].scores),[['value',source.columns.scores[0][0]]]);assert.equal(JSON.stringify(source),before);assert.throws(()=>{d.poses[0].scores.value=5;});assert(!('columns' in d));
 const result=H.applyAnalysis(d,{pose_ids:['real-id-b','real-id-a'],xy:[[1,2],[3,4]],labels:[2,1],clusters:[],projection:{}});assert.equal(result.poses[0].scores.value,d.poses[0].scores.value);assert.deepEqual(result.poses[0].contacts,d.poses[0].contacts);assert.equal(typeof Object.getOwnPropertyDescriptor(result.poses[0],'contacts').get,'function');
});
test('unknown versions, malformed dimensions/indices, mixed revisions and identities reject',async()=>{
 const changes=[d=>d.columns.version=2,d=>d.format='columns-v2',d=>d.columns.snapshot='other',d=>d.columns.state='raw',d=>d.columns.pose_ids.reverse(),d=>d.columns.dimensions.poses=3,d=>d.columns.scores[0].pop(),d=>d.columns.contacts[0].push(2),d=>d.columns.presence.scores[1][0]=1,d=>d.columns.map_presence.scores.pop(),d=>d.columns.score_keys[1]='value'];
 for(const change of changes){const d=wire();change(d);await assert.rejects(()=>T.adapt(d));}
});
test('compact rows reject mixed legacy maps regardless of presence mask',async()=>{
 for(const field of ['scores','percentiles','contacts'])for(const present of [false,true]){const d=wire();d.columns.map_presence[field][0]=present;d.poses[0][field]=field==='contacts'?['injected']:{injected:999};await assert.rejects(()=>T.adapt(d),/mixed legacy/);}
});
test('existing real rows decode exactly, including all scores and export serialization',async()=>{
 const run=spawnSync((process.env.PYTHON || 'python3'),[path.join(__dirname,'transport_fixture.py')],{encoding:'utf8',env:{...process.env,OMP_NUM_THREADS:'1',OPENBLAS_NUM_THREADS:'1',MKL_NUM_THREADS:'1',NUMEXPR_NUM_THREADS:'1'},timeout:60000,maxBuffer:4*1024*1024});assert.equal(run.status,0,run.stderr);const fixture=JSON.parse(run.stdout),d=await T.adapt(fixture.wire);assert.deepEqual(JSON.parse(JSON.stringify(d)),fixture.source);assert(Object.keys(d.poses[0].scores).length>100);assert.equal(H.csv(d.poses.map(p=>({...p.scores}))),H.csv(fixture.source.poses.map(p=>({...p.scores}))));
});
