'use strict';
const test=require('node:test'),assert=require('node:assert/strict'),path=require('node:path');
const H=require(process.env.CYP_REFERENCE_HELPERS||path.resolve(__dirname,'../static/helpers.js'));
const point={status:'AVAILABLE',xy:[1.23456789012345,-2.34567890123456]};
function fixture(){return{snapshot:'runtime',source_snapshot:'scientific',ligand:'OCNT-2313300',state:'paired',poses:[{id:'pose::raw',pose_xy:[3,4]}],atoms:[],coverage:{shown:1},revealed_reference:{id:'revealed/OCNT-2313300/1/A/LIG/601',role:'revealed_reference',snapshot:'runtime',source_snapshot:'scientific',ligand:'OCNT-2313300',state:'paired',revision:'a'.repeat(64),source_pdb_sha256:'b'.repeat(64),pose:{...point,coordinate_state:'paired',metric:'rmsd',atom_indices:[]},receptor:{...point},energy:{...point}}};}
const req={ligand:'OCNT-2313300',state:'paired',metric:'rmsd',atoms:[]};
test('reference identity and finite coordinates reject mixed or malformed packets',()=>{
 const data=fixture();assert.equal(H.checkLigand(data),data);
 for(const change of [{snapshot:'old'},{source_snapshot:'old'},{ligand:'OCNT-2312554'},{state:'raw'},{role:'pose'},{revision:'bad'},{source_pdb_sha256:null}])assert.throws(()=>H.checkLigand({...data,revealed_reference:{...data.revealed_reference,...change}}));
 for(const field of ['pose','receptor','energy'])for(const change of [{status:'AVAILABLE',xy:null},{status:'AVAILABLE',xy:[NaN,1]},{status:'UNAVAILABLE',xy:[0,0]},{status:'BAD',xy:null}])assert.throws(()=>H.checkLigand({...data,revealed_reference:{...data.revealed_reference,[field]:{...data.revealed_reference[field],...change}}}));
});
test('metric or selected atom changes cannot reuse a stale reference projection',()=>{
 const data=fixture();assert.deepEqual(H.referencePoint(data,req,'posePlot').xy,point.xy);
 for(const metric of ['contacts','heme','selected_atoms']){const p=H.referencePoint(data,{...req,metric,atoms:[0,1]},'posePlot');assert.equal(p.status,'PENDING');assert.equal(p.xy,null);}
 data.revealed_reference.pose.metric='selected_atoms';data.revealed_reference.pose.atom_indices=[0,1];
 assert.deepEqual(H.referencePoint(data,{...req,metric:'selected_atoms',atoms:[0,1]},'posePlot').xy,point.xy);
 assert.equal(H.referencePoint(data,{...req,metric:'selected_atoms',atoms:[0,2]},'posePlot').xy,null);
 assert.equal(H.referencePoint(data,{...req,state:'raw'},'energyPlot'),null);
 assert.equal(H.referencePoint(data,{...req,ligand:'OCNT-2395503'},'receptorPlot'),null);
});
test('analysis replaces the reference without modifying or adding model rows',()=>{
 const data=fixture(),before=JSON.stringify(data);
 const result={snapshot:'runtime',source_snapshot:'scientific',pose_ids:['pose::raw'],xy:[[5,6]],labels:[2],clusters:[],revealed_reference:{...data.revealed_reference,pose:{...data.revealed_reference.pose,metric:'contacts',status:'UNAVAILABLE',xy:null,reason:'not qualified'}}};
 const applied=H.applyAnalysis(data,result);assert.equal(JSON.stringify(data),before);assert.equal(applied.poses.length,1);assert.equal(applied.revealed_reference.pose.xy,null);assert.deepEqual(applied.poses[0].pose_xy,[5,6]);
 assert.equal(H.applyAnalysis(data,{...result,revealed_reference:undefined}).revealed_reference,null);
 assert.throws(()=>H.applyAnalysis(data,{...result,revealed_reference:{...result.revealed_reference,state:'raw'}}));
});
