# Explicit live API checks

These former `tests/` demonstrations contact Nebius and may incur API charges.
They are outside automated test discovery, and importing a module no longer
executes its API call. From the project root, in `biohack-agent`:

```sh
python scripts/smoke/list_nebius_models.py
python scripts/smoke/check_nebius_model.py
python scripts/smoke/check_tool_calling.py
```

Each reads the root `.env`. For the current SDK agent, prefer `python -m agent`;
the tool-calling script is the original direct-API demonstration, not the SDK
implementation. Offline tests are under `tests/agent/` and `tests/research/`.
