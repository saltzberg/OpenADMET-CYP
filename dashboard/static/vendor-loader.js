/* Classic-script loading for the vendored UMD bundle. No module exports required. */
(function(root){'use strict';
const url='/articles/admet-cyp-structure/cofolding/static/vendor/molstar.js';
const maxAttempts=3, timeoutMs=15000;
let flight=null, failure=null, attempts=0, state='idle';
function show(text,retry){
 const status=document.getElementById('viewerStatus');if(!status)return;
 status.replaceChildren(document.createTextNode(text));
 if(retry){const button=document.createElement('button');button.type='button';button.textContent='Retry viewer';button.style.pointerEvents='auto';button.onclick=()=>warm({retry:true}).catch(()=>{});status.append(' ',button);}
}
function warm({retry=false}={}){
 if(root.molstar?.Viewer){if(state!=='ready')flight=Promise.resolve(root.molstar);state='ready';return flight||(flight=Promise.resolve(root.molstar));}
 if(state==='loading')return flight;
 if(failure&&!retry)return flight;
 if(attempts>=maxAttempts)return Promise.reject(failure||new Error('Viewer retry limit reached; reload the page.'));
 attempts++;state='loading';failure=null;show('Loading viewer…',false);
 flight=new Promise((resolve,reject)=>{
  const script=document.createElement('script');script.src=url;script.async=true;
  let settled=false;
  const finish=error=>{if(settled)return;settled=true;clearTimeout(timer);script.onload=script.onerror=null;
   if(error){script.remove();state='failed';failure=error;show('Viewer unavailable · '+error.message,attempts<maxAttempts);reject(error);}
   else{state='ready';show('',false);resolve(root.molstar);}
  };
  const timer=setTimeout(()=>finish(new Error('Viewer bundle timed out')),timeoutMs);
  script.onload=()=>finish(root.molstar?.Viewer?null:new Error('Viewer bundle has no molstar.Viewer global'));
  script.onerror=()=>finish(new Error('Viewer bundle could not load'));
  document.head.append(script);
 });
 return flight;
}
root.ViewerVendor=Object.freeze({warm,get state(){return state;},get attempts(){return attempts;},maxAttempts,timeoutMs});
// Plotly stays in the ordered defer list. Consumers can use this explicit check.
root.PlotlyVendor=Object.freeze({ready:()=>root.Plotly?Promise.resolve(root.Plotly):Promise.reject(new Error('Plotly bundle unavailable'))});
})(globalThis);
