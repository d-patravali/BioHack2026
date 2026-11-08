# Biomedical research agent

A single OpenAI Agents SDK agent using **Nebius Token Factory**, with one local
`get_subtypes(disease)` tool. No KG, GNN, prediction service, or frontend is connected.

From the project root:

```sh
conda activate biohack-agent
python -m pip install -r requirements.txt
python -m agent --check
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
turns do not replace successful conversation history. Runs have a six-turn
limit; the client uses a 60-second timeout and at most one retry.

## Verify routing and test

```sh
python -m unittest discover -s tests -p 'test_agent_framework.py' -v
```

These offline tests use a mock HTTP transport to execute the real SDK model →
tool → model loop. They assert every request goes to
`https://api.tokenfactory.nebius.com/v1/chat/completions` with the configured
model, and inspect the tool result returned to the model. They also cover
configuration, unsupported diseases, malformed arguments, and debug redaction.
Run this exact test pattern: older `tests/test_nebius.py`, `test_model.py`, and
`test_tool_calling.py` are live API scripts that make requests at import time.

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
with the binary `cryptography==46.0.3` wheel, and the existing root dependency
pin was updated to reproduce that working environment. Other dependency pins
and existing tests were preserved.

To extend later, add decorated functions in `tools.py` (or a dedicated module)
and register them in `create_agent`. Keep graph and GNN services behind those
tools. The current subtype catalog is static, not a live graph query.
