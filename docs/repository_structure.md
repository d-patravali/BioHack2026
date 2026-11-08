# Repository organization and migration

The reorganization separates application code, research code, notebooks, raw
data, derived data, experiment outputs, and written findings. `agent/` stays at
the root so `python -m agent` and root `.env` resolution continue to work.
`src/primekg/` remains the reusable preprocessing package. Existing research
algorithms, numerical outputs, and notebook outputs were not regenerated.

## Old → new paths

| Previous location | Current location |
| --- | --- |
| `kg.csv` | `data/raw/primekg/kg.csv` |
| `data/processed/*.csv` | `data/processed/primekg/*.csv` |
| Root `PrimeKG_ECoGraph_merge.ipynb` | `notebooks/integration/PrimeKG_ECoGraph_merge.ipynb` |
| Root `PrimeKG_ECoGraph_GNN.ipynb` | `notebooks/models/colorectal/PrimeKG_ECoGraph_GNN.ipynb` |
| Root `PrimeKG_HER2_Breast_Cancer_GNN.ipynb` | `notebooks/models/breast/PrimeKG_HER2_Breast_Cancer_GNN.ipynb` |
| Root `primekg_findings.md` | `docs/research/primekg_findings.md` |
| Previous research-focused root `README.md` | `docs/research/primekg_exploration.md` |
| Root profiling files in `reports/` | `reports/primekg/` |
| `analysis/results/` | `reports/compatibility/` |
| `reports/inference/` | Same location |
| Root agent dependency pins | `requirements/agent.txt` |
| `tests/test_agent_framework.py` | `tests/agent/test_agent_framework.py` |
| `tests/test_inference_metrics.py` | `tests/research/test_inference_metrics.py` |
| `analysis/test_compatibility.py` | `tests/research/test_compatibility.py` |
| `tests/test_nebius.py` | `scripts/smoke/list_nebius_models.py` |
| `tests/test_model.py` | `scripts/smoke/check_nebius_model.py` |
| `tests/test_tool_calling.py` | `scripts/smoke/check_tool_calling.py` |
| Notebook `kg_data/primekg_ecograph/` assumption | `data/processed/primekg_ecograph/` |
| Notebook CRC prediction output in graph folder | `artifacts/predictions/colorectal/` |
| Notebook `kg_data/her2_breast_cancer/` assumption | Raw GDSC in `data/external/gdsc/`; results in `artifacts/predictions/breast/` |

The notebook filenames on disk use **PrimeKG**, and those actual names were
preserved. Notebook input/output path configuration was adapted to the new
layout, including the existing local ECoGraph matrix and Neptune exports.
The notebooks locate the project root from the working directory or an ancestor.

## Upstream ECoGraph checkout

`ECoGraph/` is a separate Git repository recorded in the parent index as a Git
link. It intentionally stays in place; no upstream files or nested Git metadata
were modified. This checkout currently has no root `.gitmodules`, so do not
assume a fresh parent clone will automatically fetch ECoGraph. Supply the
upstream checkout separately at `ECoGraph/` using its own repository instructions.
Formalizing submodule management is a separate Git decision.

## Data and result preservation

The PrimeKG raw CSV was moved within the filesystem, not copied or rewritten.
Its checksum and all relocated derived CSV/JSON/TXT artifacts were checked
against pre-move snapshots. `ECoGraph/` and the root `.env` stayed untouched.
Historical source manifests, fingerprints, and CSV fields may still name the
old paths. These are records of prior runs and have deliberately not been
rewritten. Use the mapping above to locate the same inputs and results.

Readable report links and reproduction commands were updated. New script runs
default to the current paths, and naturally record those paths in new provenance.
Re-running analysis can overwrite its results; it is unnecessary for viewing
the existing findings.

## Component boundaries

- `agent/`: provider access, conversation loop, and registered tool interfaces.
- `src/primekg/` and `scripts/`: initial preprocessing and reproducible baselines.
- `analysis/`: compatibility investigation; outputs belong under `reports/`.
- `notebooks/`: teammate graph/model experiments awaiting scientific/code review.
- `data/processed/`: graph tables; `artifacts/`: predictions and model files.
- `reports/`: saved analysis outputs; `docs/research/`: narrative research context.

This layout does not register notebook outputs as agent tools, validate the GNN,
change biological relationship semantics, or choose a final modeling architecture.
