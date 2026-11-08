# Hackathon demonstration pack

Start with [pitch_outline.md](pitch_outline.md) and the
[interactive prompts](../../agent/DEMO.md).

Both prompts were rehearsed in one live Nebius conversation. The backup is a
recording, not fresh inference:

- [Recorded walkthrough](recordings/20261108_042554/walkthrough.md)
- [Structured tool calls and outputs](recordings/20261108_042554/transcript.json)

Offline replay, with no API calls or credentials:

```sh
python3 scripts/rehearse_agent.py --replay docs/demo/recordings/20261108_042554/transcript.json
```

Start the interactive agent:

```sh
conda activate biohack-agent
python -m agent --check-data
python -m agent
```

For another recorded live run, use `python scripts/rehearse_agent.py --live` in
the agent environment. Each run creates a new directory; existing recordings
are never overwritten. Live runs can incur Nebius charges. The recorder keeps
tool result payloads and model text, redacts the configured key, and fingerprints
six graph/experiment artifacts before and after. No notebooks or models run.

## Review notes

Tool-evidence checks passed for both scenarios. The source data/results were
unchanged. This is a smoke test, not a systematic test of agent truthfulness.

The second recorded answer's wording, “neither APC's profile counts nor its HPO
associations explain the tied drug ranking,” should be presented carefully:
**mutation frequencies are model inputs; HPO annotations are not.** Neither
alone establishes a causal/faithful attribution for the ranking. Agent
instructions were clarified after this rehearsal. The transcript retains the
actual recorded response rather than rewriting it to match the correction.

The results do not demonstrate GNN superiority, clinical efficacy, literature
novelty, reduced hallucination rates, or a phenotype-aware trained GNN. Present
the engineering prototype and its measured limitations directly.
