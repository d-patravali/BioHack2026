"use strict";
const state={mode:"recorded",backend:"ecograph",ready:false,busy:false,session:crypto.randomUUID(),info:null,traces:new Map()};
const $=id=>document.getElementById(id);
const esc=s=>String(s??"").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const inline=s=>esc(s).replace(/\*\*([^*]+)\*\*/g,"<strong>$1</strong>").replace(/`([^`]+)`/g,"<code>$1</code>");
function markdown(text){let result="",inList=false;for(const line of String(text).split("\n")){if(/^\s*[-*]\s+/.test(line)||/^\s*\d+[.)]\s+/.test(line)){if(!inList){result+="<ul>";inList=true}result+="<li>"+inline(line.replace(/^\s*(?:[-*]|\d+[.)])\s+/,""))+"</li>"}else{if(inList){result+="</ul>";inList=false}if(/^#+\s/.test(line))result+="<h3>"+inline(line.replace(/^#+\s*/,""))+"</h3>";else if(line.trim())result+="<p>"+inline(line)+"</p>"}}return result+(inList?"</ul>":"")}
function banner(message){$("banner").hidden=!message;$("banner").textContent=message||""}
function busy(value){
  state.busy=value;
  $("send").disabled=value||state.mode!=="live"||!state.ready;
  $("prompt").disabled=$("send").disabled;
  $("reset").disabled=value;
  $("backend").disabled=value||state.mode!=="live"||!state.ready;
  document.querySelectorAll("[data-scenario],.mode-switch button").forEach(button=>{
    button.disabled=value||!state.ready||(button.id==="mode-live"&&!state.info?.live_available);
  });
  $("connection-status").textContent=value?"Loading…":state.mode==="live"?`${state.backend==="crcsc"?"CRCSC":"ECoGraph"} · Nebius ready`:"ECoGraph · offline replay ready";
}

function scroll(){const pane=$("messages");pane.scrollTop=pane.scrollHeight}
function addMessage(role,text){$("empty-state")?.remove();const box=document.createElement("article");box.className="message "+role;box.innerHTML=`<div class="message-meta">${role==="user"?"YOU":state.mode==="live"?`RESEARCH AGENT · LIVE · ${state.backend==="crcsc"?"CRCSC":"ECoGraph"}`:"RESEARCH AGENT · RECORDED · ECoGraph"}</div><div class="message-body"></div>`;box.querySelector(".message-body").innerHTML=markdown(text);$("messages").append(box);scroll();return box}
function summary(name,result){if(result?.status!=="ok")return result?.message||result?.status||"Tool returned";if(name==="get_drug_predictions")return `${result.returned_count} candidates · full cutoff ties preserved · cached scores`;if(name==="get_drug_targets")return `${result.target_count} recorded direct target(s) for ${result.drug_name}`;if(name==="get_subtype_genes")return `${result.returned_count} of ${result.total_genes} ${result.profile_type==="expression_signature"?"expression signature":"retained mutation-profile"} genes`;if(name==="get_gene_phenotypes")return `${result.phenotype_count} recorded direct associations for NCBI ${result.gene_id}`;if(name==="get_model_evaluation")return "Popularity and molecular-context ablation retrieved";return "Source returned successfully"}
function handle(item,answerBox){if(item.type==="tool_call"){const trace=document.createElement("div");trace.className="trace running";trace.innerHTML=`<div class="trace-head"><strong>${esc(item.tool)}</strong><span>Querying</span></div><div class="trace-summary">${esc(JSON.stringify(item.arguments))}</div><details><summary>Inspect tool data</summary><pre></pre></details>`;answerBox.before(trace);state.traces.set(item.call_id,trace);scroll()}else if(item.type==="tool_result"){const trace=state.traces.get(item.call_id);if(trace){trace.classList.remove("running");trace.querySelector(".trace-head span").textContent=item.result?.status==="ok"?"Complete":"Check result";trace.querySelector(".trace-summary").textContent=summary(item.tool,item.result);trace.querySelector("pre").textContent=JSON.stringify(item.result,null,2)}scroll()}else if(item.type==="answer"){answerBox.querySelector(".message-body").innerHTML=markdown(item.text);if(item.review_note){const note=document.createElement("div");note.className="review-note";note.textContent="Presenter note: "+item.review_note;answerBox.append(note)}scroll()}else if(item.type==="error"){answerBox.querySelector(".message-body").textContent=item.message;banner(item.message)}else if(item.type==="start"&&item.mode==="recorded"){const date=new Date(item.recorded_utc).toLocaleString();$("connection-status").textContent="Recorded "+date}}
async function loadBackend(backend){
  const response=await fetch(`/api/bootstrap?backend=${encodeURIComponent(backend)}`);
  if(!response.ok)throw new Error("Could not load the selected backend");
  const info=await response.json();
  state.info=info;state.backend=backend;state.ready=true;
  $("backend").value=backend;
  $("subtype-label").textContent=info.selected_backend.subtypes.join(" · ");
  $("gnn-label").textContent=info.selected_backend.gnn_label;
  const conditions=info.evaluation.conditions;
  $("popularity").textContent=conditions?.drug_popularity?.mean_test_ap?.toFixed(4)||"—";
  $("gnn").textContent=conditions?.[info.selected_backend.gnn_condition]?.mean_test_ap?.toFixed(4)||"—";
  $("ablation").textContent=conditions?.drug_only?.mean_test_ap?.toFixed(4)||"—";
  if(!conditions)banner("Selected backend results are unavailable. No alternative experiment has been substituted.");
  for(const scenario of info.scenarios){
    document.querySelectorAll(`[data-scenario="${scenario.id}"]`).forEach(button=>{
      const title=button.querySelector(".scenario-title"),description=button.querySelector("small");
      if(title)title.textContent=scenario.short_title;
      if(description)description.textContent=scenario.description;
    });
  }
  $("provider").textContent=info.model==="Not configured"?"Live agent not configured · ECoGraph replay available":`Nebius · ${info.model}`;
}
function renderMode(){
  const live=state.mode==="live";
  $("mode-recorded").classList.toggle("selected",!live);
  $("mode-live").classList.toggle("selected",live);
  $("mode-recorded").setAttribute("aria-pressed",String(!live));
  $("mode-live").setAttribute("aria-pressed",String(live));
  $("mode-badge").textContent=live?`LIVE · ${state.backend.toUpperCase()}`:"RECORDED · ECOGRAPH";
  $("mode-badge").classList.toggle("live",live);
  $("mode-description").textContent=live?"Live conversation; graph and model results stay cached.":"Replay the original ECoGraph session offline.";
  $("backend-hint").textContent=live?"Changing backend starts a new conversation.":"Recordings use ECoGraph. Switch to Live agent to choose CRCSC.";
  $("composer-mode").textContent=live?`${state.info?.selected_backend.label} · cached GNN`:"ECoGraph recording · no fresh inference";
  $("prompt").placeholder=live?(state.backend==="crcsc"?"Ask about CMS drug hypotheses or expression signatures…":"Ask about GI subtype drugs or mutation observations…"):"Choose a recorded investigation, or switch to Live agent…";
}
async function reset(force=false){
  if(state.busy&&!force)return;
  try{await fetch("/api/reset",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({session_id:state.session})})}catch{}
  state.session=crypto.randomUUID();
  $("messages").innerHTML='<div id="empty-state" class="empty-state"><span class="empty-symbol">⌘</span><h2>A ranking is the start of a question.</h2><p>Choose a guided investigation, or ask a question in Live mode.</p></div>';
  state.traces.clear();banner("");$("prompt").value="";
}
async function setMode(mode){
  if(state.busy||mode===state.mode)return;
  busy(true);
  try{
    await reset(true);
    if(mode==="recorded"&&state.backend!=="ecograph")await loadBackend("ecograph");
    state.mode=mode;renderMode();
  }catch(error){banner(error.message)}finally{busy(false)}
}
async function setBackend(backend){
  if(state.busy||state.mode!=="live")return;
  busy(true);
  try{await reset(true);await loadBackend(backend);renderMode()}
  catch(error){$("backend").value=state.backend;banner(error.message)}finally{busy(false)}
}
async function send(prompt,scenario){if(state.busy)return;const record=state.info?.scenarios.find(s=>s.id===scenario);prompt=record?.prompt||prompt;if(!prompt?.trim())return;banner("");state.traces.clear();addMessage("user",prompt);const answerBox=addMessage("assistant","Investigating the available evidence…");busy(true);$("prompt").value="";try{const response=await fetch("/api/chat",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({session_id:state.session,mode:state.mode,backend:state.backend,prompt,scenario_id:scenario})});if(!response.ok){const data=await response.json();throw new Error(data.error||"The request failed")}const reader=response.body.getReader(),decoder=new TextDecoder();let buffer="";while(true){const {value,done}=await reader.read();buffer+=decoder.decode(value||new Uint8Array(),{stream:!done});let end;while((end=buffer.indexOf("\n\n"))!==-1){const block=buffer.slice(0,end);buffer=buffer.slice(end+2);const line=block.split("\n").find(l=>l.startsWith("data: "));if(line)handle(JSON.parse(line.slice(6)),answerBox)}if(done)break}}catch(error){answerBox.querySelector(".message-body").textContent=error.message;banner(error.message+" You can switch to Recorded mode for the verified offline walkthrough.")}finally{busy(false);scroll()}}
document.querySelectorAll("[data-scenario]").forEach(button=>button.addEventListener("click",()=>send("",button.dataset.scenario)));
$("mode-recorded").addEventListener("click",()=>setMode("recorded"));$("mode-live").addEventListener("click",()=>setMode("live"));$("reset").addEventListener("click",()=>reset());$("backend").addEventListener("change",()=>setBackend($("backend").value));$("chat-form").addEventListener("submit",e=>{e.preventDefault();send($("prompt").value)});$("prompt").addEventListener("keydown",e=>{if(e.key==="Enter"&&!e.shiftKey){e.preventDefault();send($("prompt").value)}});
(async()=>{
  busy(true);
  try{await loadBackend("ecograph");renderMode();if(!state.info.recorded_available)banner("ECoGraph recorded files are unavailable.")}
  catch(error){banner(error.message);state.ready=false}
  finally{busy(false)}
})();
