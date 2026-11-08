# Team Rocket — local research-agent demo

A single-page interface around the existing Python agent. It runs locally; credentials stay in Python. No graph merge or GNN training occurs.

## Launch

From the repository root:

```sh
conda activate biohack-agent
python -m ui
```

Open **http://127.0.0.1:8765**. Use `python -m ui --port 8766` if the port is occupied. Stop with Ctrl+C. Dependencies are already listed in `requirements/agent.txt`; a fresh environment can install them using `python -m pip install -r requirements/agent.txt`.

## Demonstration

1. Start in **Recorded session** and choose **Investigate a 29-drug tie**. Follow the tool calls and expand their results to inspect the popularity baseline and ablation.
2. Choose **Explore APC evidence** to show subtype mutation observations and recorded phenotype–gene annotations.
3. Switch to **Live agent** to ask your own questions or rerun a guided prompt. Live mode requires `NEBIUS_API_KEY` and `NEBIUS_MODEL` in the root `.env`. It uses the existing Nebius Chat Completions configuration, never an OpenAI endpoint. Model answers can vary.

Recorded mode replays the actual saved Nebius conversation at `docs/demo/recordings/20261108_042554/transcript.json`; it makes no provider calls and accepts only the two recorded scenarios. The presenter note clarifies model inputs versus retrieved phenotype evidence. A recorded response is not a fresh model run.

Live mode streams tool start/result events and the final response. Conversation history stays in server memory until reset or shutdown. Expandable tool traces show each call and its returned data within the conversation. Switching modes resets conversation history.

## Scientific boundaries

- GNN scores come from the reviewed cached experiment, including cutoff ties. There is no fresh GNN inference.
- Raw logits are not efficacy probabilities. Drug-target connectivity is not therapeutic benefit.
- The popularity baseline outperformed the molecular-context GNN on the recorded sampled indication-recovery evaluation. Subtype efficacy has not been evaluated.
- Mutation frequencies are model inputs; phenotype edges are separately retrieved annotations, not GNN attributions.
- Research files and datasets are read only. The UI does not write transcripts or alter model outputs.

The server binds only to `127.0.0.1`. It is an unauthenticated local demo, not a hosted multi-user application; keep it local. API keys are never sent to the browser. Missing configuration disables Live mode; request failures show a safe message and suggest recorded replay.

## Verify

```sh
python -m unittest discover -s tests/agent -p 'test_*.py' -v
node --check ui/static/app.js
```

HTTP tests exercise actual recorded replay without a provider, live history/reset with an injected backend, and request validation. The SDK tests separately check that calls use Nebius Chat Completions. Tests make no paid API calls.

Files: `server.py` provides HTTP/event streaming; `__main__.py` launches Uvicorn; `static/` contains the buildless HTML/CSS/JavaScript interface. The existing `agent/` package remains the tool/model backend.

## Select a data backend

In **Live agent** mode, the **Data backend** selector chooses the original
ECoGraph GI workflow or CRCSC CMS expression workflow. Switching clears history
and updates the evaluation panel and guided prompts. Recorded mode returns to
ECoGraph because only those original recordings exist. CRCSC exposes cached
seed42 PPI/pathway signed-profile scores and expression signatures; there is no
fresh inference or literature verification. See [backend guide](../agent/BACKENDS.md).
