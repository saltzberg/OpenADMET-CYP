const assert=require('node:assert/strict');
const {createHash,randomBytes,webcrypto}=require('node:crypto');
const {fallback,digest}=require('../static/sha256.js');
(async()=>{
 for(const size of [0,1,3,55,56,63,64,65,127,128,1000,1000000]){
  const bytes=randomBytes(size),expected=createHash('sha256').update(bytes).digest('hex');
  assert.equal(Buffer.from(fallback(bytes)).toString('hex'),expected);
  assert.equal(Buffer.from(await digest('SHA-256',bytes)).toString('hex'),expected);
  const native=await webcrypto.subtle.digest('SHA-256',bytes);
  assert.equal(Buffer.from(fallback(bytes)).toString('hex'),Buffer.from(native).toString('hex'));
 }
 const padded=Buffer.from('xabcx');assert.equal(Buffer.from(fallback(padded.subarray(1,4))).toString('hex'),createHash('sha256').update('abc').digest('hex'));
 await assert.rejects(digest('SHA-1',new Uint8Array()),/Only SHA-256/);
 console.log('PASS: SHA-256 fallback matches Node and Web Crypto for 12 sizes, sliced views and unsupported-algorithm rejection.');
})().catch(e=>{console.error(e);process.exitCode=1;});
