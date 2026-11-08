# BioHack 2026 — biomedical graph research agent

This repository contains a working Nebius-backed agent framework, dataset
compatibility research, and teammate notebooks for graph integration and GNN
experiments. The agent exposes a subtype catalog, reviewed cached GNN outputs,
local graph evidence, and model evaluation/ablation results. It does not perform
fresh GNN inference. See the [research-agent walkthrough](agent/DEMO.md).

## Find the right component

```text
agent/                         Reusable CLI agent, configuration, and tools
src/primekg/                   Reusable PrimeKG preprocessing
analysis/                      Compatibility analysis code
scripts/                       Preprocessing and baseline entry points
  smoke/                       Explicit live Nebius API checks
notebooks/
  integration/                 PrimeKG–ECoGraph merge notebook
  models/colorectal/            CRC indication/subtype GNN experiment
  models/breast/                HER2 breast-cancer response experiment
data/
  raw/primekg/kg.csv            Original PrimeKG CSV (ignored by Git)
  external/                    Downloaded HGNC/GDSC references (ignored)
  processed/                   Derived tables and future merged graphs (ignored)
ECoGraph/                      Separate upstream repository; kept intact
artifacts/                     Future model predictions and checkpoints (ignored)
reports/
  primekg/                     Dataset profiling results
  inference/                   Initial annotation-recovery baselines
  compatibility/               ECoGraph–PrimeKG compatibility results
docs/                          Repository guide and research writeups
requirements/                  Separate agent, research, notebook dependencies
tests/
  agent/                       Offline SDK/tool tests
  research/                    Offline analysis/metric tests
```

See [repository structure and migration notes](docs/repository_structure.md)
for exact old-to-new paths and component boundaries.

## Agent

```sh
conda activate biohack-agent
python -m pip install -r requirements/agent.txt
python -m agent --check
python -m agent
```

Keep `NEBIUS_API_KEY` and `NEBIUS_MODEL` in the root `.env`, using
[.env.example](.env.example) as a template. See [agent instructions](agent/README.md).
`--check` is offline. Actual prompts contact Nebius; OpenAI tracing is disabled.

## Research workflows

Run commands from the repository root. Python defaults resolve relocated data
and output directories relative to the checkout; explicit CLI overrides remain available.

```sh
python3 -m pip install -r requirements/research.txt
python3 scripts/prepare_primekg.py
python3 scripts/explore_inference.py
python3 analysis/run_compatibility.py
```

Existing results can be read without rerunning large scans:

- [PrimeKG exploration](docs/research/primekg_exploration.md)
- [Biological inference findings](docs/research/primekg_findings.md)
- [Compatibility report](reports/compatibility/compatibility_report.md)
- [Compatibility commands and methods](analysis/README.md)
- [Incremental CRC model experiments](analysis/CRC_EXPERIMENTS.md): weighted
  subtype profiles, molecular ablations, baseline comparisons and tie audits.

Stored provenance files describe the original execution layout; consult the
migration notes to translate historical paths. Their hashes and numerical
results were preserved rather than regenerated during organization.

## Modeling notebooks

See [notebook workflow](notebooks/README.md). The intended dependency order is:

1. [PrimeKG–ECoGraph integration](notebooks/integration/PrimeKG_ECoGraph_merge.ipynb)
2. [CRC GNN experiment](notebooks/models/colorectal/PrimeKG_ECoGraph_GNN.ipynb)
3. [Breast-cancer experiment](notebooks/models/breast/PrimeKG_HER2_Breast_Cancer_GNN.ipynb),
   which also reads the exported graph tables.

The notebooks are research code awaiting review, not deployed services.
Relocation updated paths only; no notebook was executed or its model validated.

## Dependencies and offline tests

`requirements.txt` includes the agent and lightweight research requirements.
Notebook/GNN packages are optional in `requirements/notebooks.txt`; install them
in a modeling environment rather than adding PyTorch to the minimal agent environment.

```sh
# Agent environment:
python -m unittest discover -s tests/agent -p 'test_*.py' -v
# Research environment:
python3 -m unittest discover -s tests/research -p 'test_*.py' -v
# An environment with all lightweight requirements and pytest:
python -m pytest
```

Live checks live under `scripts/smoke/` and run only when explicitly invoked.
Pytest discovery is restricted to `tests/`. No test requires the merged graph,
a trained GNN, or live credentials.
