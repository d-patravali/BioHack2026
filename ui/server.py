"""Buildless local UI; all credentials and tool execution stay in Python."""

import asyncio
import json
import re
from pathlib import Path

from starlette.applications import Starlette
from starlette.middleware import Middleware
from starlette.middleware.trustedhost import TrustedHostMiddleware
from starlette.requests import Request
from starlette.responses import FileResponse, JSONResponse, StreamingResponse
from starlette.routing import Route

from agent.config import PROJECT_ROOT, load_settings
from agent.repository import query
from agent.backends import BACKENDS, catalog, make_repository, validate_backend
from agent.runtime import DebugHooks, create_agent, create_client, run_turn
from agent.__main__ import error_message
from scripts.rehearse_agent import SCENARIOS

STATIC = Path(__file__).resolve().parent / "static"
RECORDING = PROJECT_ROOT / "docs/demo/recordings/20261108_042554/transcript.json"


def configuration():
    try:
        settings = load_settings()
        return {"live_available": True, "model": settings.model.replace(settings.api_key, "[REDACTED]")}
    except ValueError:
        return {"live_available": False, "model": "Not configured"}


def event(payload):
    return "data: " + json.dumps(payload, ensure_ascii=False) + "\n\n"


def read_recording():
    try:
        return json.loads(RECORDING.read_text())
    except (OSError, ValueError):
        return None


async def live_answer(prompt, history, emit, backend='ecograph'):
    settings = load_settings()
    class BrowserHooks(DebugHooks):
        def cleaned(self, value):
            try:
                return json.loads(self.safe(value))
            except (ValueError, TypeError):
                return self.safe(value)
        async def on_tool_start(self, context, agent, tool):
            await emit({"type": "tool_call", "tool": tool.name,
                "call_id": getattr(context, "tool_call_id", tool.name),
                "arguments": self.cleaned(getattr(context, "tool_arguments", "{}"))})
        async def on_tool_end(self, context, agent, tool, result):
            await emit({"type": "tool_result", "tool": tool.name,
                "call_id": getattr(context, "tool_call_id", tool.name), "result": self.cleaned(result)})
    hooks = BrowserHooks(settings.api_key)
    async with create_client(settings) as client:
        result = await run_turn(create_agent(settings, client, backend=backend), history, prompt, hooks)
    return hooks.safe(result.final_output), result.to_input_list()


def scenarios(backend):
    if backend == 'crcsc':
        return [
            {'id': 'candidate_audit', 'title': 'Inspect CMS2 drug hypotheses',
             'short_title': 'Why do these drugs tie?', 'description': 'Inspect CMS2 targets and model controls',
             'prompt': 'Retrieve the top five cached CMS2 drug hypotheses, inspect their shared targets and model evaluation, and explain the ties and limitations. Do not claim efficacy or literature verification.'},
            {'id': 'phenotype_gene_evidence', 'title': 'Inspect CMS1 expression context',
             'short_title': 'Inspect CMS expression', 'description': 'Separate expression from mutation evidence',
             'prompt': 'Retrieve the three largest-weight CMS1 expression signature genes. Report higher/lower expression, signed weights and discovery FDR; retrieve phenotype context for the first gene. Separate recorded context from causal or therapeutic evidence.'}]
    titles = [('Why do these drugs tie?', 'Inspect candidates, targets & controls'),
              ('Trace phenotype evidence', 'Connect APC observations & annotations')]
    return [{**s, 'short_title': titles[i][0], 'description': titles[i][1]} for i,s in enumerate(SCENARIOS)]


def create_app(answer=live_answer):
    repositories = {name: make_repository(name) for name in BACKENDS}
    sessions = {}

    async def home(request):
        return FileResponse(STATIC / "index.html")

    async def asset(request):
        name = request.path_params["name"]
        if name not in {"app.js", "style.css"}:
            return JSONResponse({"error": "Not found"}, status_code=404)
        return FileResponse(STATIC / name)

    async def bootstrap(request):
        backend = request.query_params.get('backend', 'ecograph')
        try:
            validate_backend(backend)
        except ValueError:
            return JSONResponse({'error': 'Choose ecograph or crcsc'}, status_code=400)
        recording = read_recording() if backend == 'ecograph' else None
        evaluation = await asyncio.to_thread(query, repositories[backend], 'get_model_evaluation')
        return JSONResponse({**configuration(), 'provider': 'Nebius Token Factory',
            'backend': backend, 'backend_catalog': [catalog(b) for b in BACKENDS],
            'selected_backend': catalog(backend), 'recorded_available': recording is not None,
            'recorded_utc': recording.get('recorded_utc') if recording else None,
            'evaluation': evaluation, 'scenarios': scenarios(backend)})

    async def health(request):
        return JSONResponse({"status": "ok", "service": "local_research_demo"})

    def validate_origin(request):
        origin = request.headers.get("origin")
        if origin and origin.rstrip("/") != str(request.base_url).rstrip("/"):
            return JSONResponse({"error": "Cross-origin requests are not accepted"}, status_code=403)
        if not request.headers.get("content-type", "").startswith("application/json"):
            return JSONResponse({"error": "Expected application/json"}, status_code=415)

    async def reset(request):
        rejection = validate_origin(request)
        if rejection:
            return rejection
        try:
            body = await request.json()
        except ValueError:
            return JSONResponse({"error": "Invalid JSON"}, status_code=400)
        if not isinstance(body, dict):
            return JSONResponse({"error": "Expected an object"}, status_code=400)
        if not isinstance(body.get("session_id"), str):
            return JSONResponse({"error": "Invalid session identifier"}, status_code=400)
        session = sessions.get(body.get("session_id"))
        if session and session["lock"].locked():
            return JSONResponse({"error": "Wait for the current turn before resetting"}, status_code=409)
        sessions.pop(body.get("session_id"), None)
        return JSONResponse({"status": "ok"})

    async def chat(request: Request):
        rejection = validate_origin(request)
        if rejection:
            return rejection
        try:
            body = await request.json()
        except ValueError:
            return JSONResponse({"error": "Invalid JSON"}, status_code=400)
        if not isinstance(body, dict):
            return JSONResponse({"error": "Expected an object"}, status_code=400)
        backend = body.get('backend', 'ecograph')
        try:
            validate_backend(backend)
        except ValueError:
            return JSONResponse({'error': 'Choose ecograph or crcsc'}, status_code=400)
        if body.get('mode') == 'recorded' and backend != 'ecograph':
            return JSONResponse({'error': 'Recorded walkthroughs belong to ECoGraph. Select Live agent for CRCSC.'}, status_code=400)
        session_id, mode = body.get("session_id"), body.get("mode")
        if not isinstance(session_id, str) or not re.fullmatch(r"[A-Za-z0-9-]{1,64}", session_id):
            return JSONResponse({"error": "Invalid session identifier"}, status_code=400)
        if not isinstance(mode, str) or mode not in {"live", "recorded"}:
            return JSONResponse({"error": "Choose live or recorded mode"}, status_code=400)
        scenario_id = body.get("scenario_id")
        prompt = body.get("prompt", "")
        scenario = next((s for s in scenarios(backend) if s["id"] == scenario_id), None)
        if scenario_id is not None and scenario is None:
            return JSONResponse({"error": "Unknown walkthrough"}, status_code=400)
        if scenario:
            prompt = scenario["prompt"]
        if not isinstance(prompt, str) or not 1 <= len(prompt.strip()) <= 4000:
            return JSONResponse({"error": "Enter a prompt of 1–4000 characters"}, status_code=400)
        if mode == "recorded":
            recording = read_recording()
            recorded = next((s for s in (recording or {}).get("scenarios", []) if s["id"] == scenario_id), None)
            if recorded is None:
                return JSONResponse({"error": "Select an available recorded walkthrough"}, status_code=400)
            async def replay():
                yield event({"type": "start", "mode": "recorded", "recorded_utc": recording["recorded_utc"]})
                for item in recorded["events"]:
                    yield event({"type": item["event"], **{k: v for k, v in item.items() if k != "event"}})
                    await asyncio.sleep(.08)
                yield event({"type": "answer", "text": recorded.get("answer", "No recorded answer available"),
                    "mode": "recorded", "review_note": "Mutation frequencies are model inputs; phenotype edges are retrieved separately. Neither alone is a faithful attribution of a drug score."})
                yield event({"type": "done"})
            return StreamingResponse(replay(), media_type="text/event-stream", headers={"Cache-Control": "no-store"})
        if not configuration()["live_available"]:
            return JSONResponse({"error": "Live mode needs NEBIUS_API_KEY and NEBIUS_MODEL in the root .env. Recorded mode is available offline."}, status_code=503)
        if session_id not in sessions:
            if len(sessions) >= 32:
                return JSONResponse({"error": "Local session limit reached. Reset an existing conversation or restart the server."}, status_code=503)
            sessions[session_id] = {"history": [], "lock": asyncio.Lock(), "backend": backend}
        session = sessions[session_id]
        if session['backend'] != backend:
            return JSONResponse({'error': 'Reset the conversation before changing data backends.'}, status_code=409)
        if session["lock"].locked():
            return JSONResponse({"error": "A turn is already running in this conversation"}, status_code=409)
        await session["lock"].acquire()

        async def stream():
            queue = asyncio.Queue()
            async def run():
                try:
                    text, history = await answer(prompt.strip(), session["history"], queue.put, backend=backend)
                    session["history"] = history
                    await queue.put({"type": "answer", "text": text, "mode": "live", "backend": backend})
                except Exception as error:
                    await queue.put({"type": "error", "message": error_message(error)})
                finally:
                    await queue.put({"type": "done"})
            task = asyncio.create_task(run())
            try:
                yield event({"type": "start", "mode": "live", "backend": backend})
                while True:
                    item = await queue.get()
                    yield event(item)
                    if item["type"] == "done":
                        break
            finally:
                if not task.done():
                    task.cancel()
                try:
                    await task
                except asyncio.CancelledError:
                    pass
                session["lock"].release()
        return StreamingResponse(stream(), media_type="text/event-stream", headers={"Cache-Control": "no-store"})

    return Starlette(routes=[Route("/", home), Route("/static/{name}", asset),
        Route("/health", health), Route("/api/bootstrap", bootstrap),
        Route("/api/chat", chat, methods=["POST"]), Route("/api/reset", reset, methods=["POST"])],
        middleware=[Middleware(TrustedHostMiddleware, allowed_hosts=["127.0.0.1", "localhost", "testserver"])],
        max_body_size=16384)


app = create_app()
