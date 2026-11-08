# Biomedical research agent

A single OpenAI Agents SDK agent using **Nebius Token Factory**, connected to
reviewed cached CRC GNN experiment outputs and read-only local graph evidence.
There is no live GNN inference service or frontend. See [DEMO.md](DEMO.md).

From the project root:

```sh
conda activate biohack-agent
python -m pip install -r requirements/agent.txt
python -m agent --check
python -m agent --check-data
python -m agent
```

The existing root `.env` must contain `NEBIUS_API_KEY` and `NEBIUS_MODEL`; see
`.env.example` for names. Never commit `.env`. Existing environment variables
override `.env` values. The root `.env` is resolved relative to this package,
not guessed from another working directory. No `OPENAI_API_KEY` is needed.

Example prompt: `What colorectal cancer subtypes are available?`

Single-prompt mode (makes a live, potentially billable Nebius request):

```sh
python -m agent --prompt "Use get_subtypes to list the colorectal cancer subtypes."
```

The model requests a function call; the SDK validates arguments, executes the
decorated local function, sends its result back to the model, and returns the
final answer. Run hooks print the tool name, arguments and result while the tool
executes. Unknown diseases return an explicit unsupported status and an empty
list. `/reset` clears in-memory conversation history; `/quit` exits. Failed
turns do not replace successful conversation history. Runs have a ten-turn
limit; the client uses a 60-second timeout and at most one retry.

## Verify routing and test

```sh
python -m unittest discover -s tests/agent -p 'test_*.py' -v
```

These offline tests use a mock HTTP transport to execute the real SDK model →
tool → model loop. They assert every request goes to
`https://api.tokenfactory.nebius.com/v1/chat/completions` with the configured
model, and inspect the tool result returned to the model. They also cover
configuration, unsupported diseases, malformed arguments, and debug redaction.
Live API demonstrations are in `scripts/smoke/` and execute only when invoked
explicitly. Automated tests stay offline.

`--check` prints the fixed provider URL and configured model and constructs the
agent offline. `runtime.py` explicitly constructs `OpenAIChatCompletionsModel`
with an `AsyncOpenAI` client whose base URL is fixed to Nebius, rather than a
default model name resolved by the SDK. OpenAI tracing is disabled globally and
per run, so the SDK does not export traces to OpenAI. OpenAI endpoint environment
overrides are not used. Keys are excluded from settings repr and redacted from
CLI/debug output; provider exception bodies are not printed. Tool debugging is
intentionally visible, so avoid entering sensitive patient information.

The implementation follows the installed SDK 0.24.0 APIs and the official
[models and providers documentation](https://developers.openai.com/api/docs/guides/agents/models).

The existing Conda installation had a `cryptography` 50.0.2/OpenSSL dynamic-link
error (`_EVP_DigestSqueeze`) that prevented importing `agents`. It was repaired
with the binary `cryptography==46.0.3` wheel, and the agent dependency
pin in `requirements/agent.txt` was updated to reproduce that working environment. Other dependency pins
and existing tests were preserved.

To extend later, add decorated functions in `tools.py` and register them in
`create_agent`. `repository.py` contains the read-only backend adapters, with
standard-library CSV/gzip readers and lazy in-process caches. No pandas, NumPy,
PyTorch, or new environment installs are required by the agent. Missing files
return `backend_unavailable` instead of fabricated results. The subtype catalog
is static; phenotype/target tools query tables derived from PrimeKG.

Tools now include `get_subtypes`, `get_drug_predictions`, `get_drug_targets`,
`get_subtype_genes`, `get_gene_phenotypes`, and `get_model_evaluation`.
Predictions pin reviewed run `20261108_034016`, target-only model seed42,
frequency-weighted profile. They are cached, not recomputed. All candidates tied
at the requested cutoff are returned; `top_k=5` currently returns a complete
29-drug tie. Scores are raw logits, not probabilities. The actual seed42 test AP
is 0.8158, distinct from the three-seed mean 0.8254; popularity AP is 0.8941.
The tool returns both selected-model and aggregate evaluation metadata.

`get_subtype_genes` supports a `gene` filter across the whole retained profile,
so absence from a displayed top list is not confused with absence from the full
profile. `get_gene_phenotypes` queries only explicit `phenotype_protein` edges,
deduplicates reciprocal edges and returns unknown/ambiguous inputs explicitly.
The trained CRC model does not use phenotype edges. Those tool results are
retrieved context, not a faithful model explanation or efficacy evidence.

`--check-data` checks all local backends without loading `.env` or contacting
Nebius. Data tests use synthetic fixtures; the SDK integration test exercises
prediction -> targets/phenotypes/subtype evidence -> answer with a mock transport.
