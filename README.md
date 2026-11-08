# Team Rocket — BioHack 2026

An agent-assisted biomedical knowledge-graph research prototype built around the
hackathon prompt:

> How can we prototype agentic frameworks that leverage GNNs to parse phenotype-genotype knowledge graphs?

The application lets researchers investigate **cached drug–subtype association
hypotheses**, inspect recorded molecular and phenotype–gene evidence, and compare
GNN results with simple baselines. It includes a local web interface, a terminal
agent, and reproducible research experiments.

**This is a research proof of concept.** The agent does not run fresh GNN
inference, establish drug efficacy, or perform literature verification. The
trained drug models did not use phenotype edges; the agent retrieves those
annotations separately.

## Quick start: local demo

Run from the repository root using the existing Python 3.11 environment:

```sh
conda activate biohack-agent
python -m pip install -r requirements/agent.txt
python -m ui
```

If the environment does not exist, first create it with
`conda create -n biohack-agent python=3.11`.

Open **http://127.0.0.1:8765**. Keep the server terminal running; stop it with
Ctrl+C. Use `python -m ui --port 8766` if the port is occupied. If the application
is already running, refresh the page rather than start a second server.

### Recorded mode — no API connection required

1. Select **Recorded**.
2. Choose **ECoGraph · GI genomic (original)** or **CRCSC · CMS expression** under
   **Data backend**.
3. Click a guided investigation:
   - **Why do these drugs tie?** inspects a saved ranking, targets, and evaluation.
   - **Trace phenotype evidence** (ECoGraph) or **Inspect CMS expression** (CRCSC)
     inspects subtype genes and recorded phenotype context.

These are actual saved Nebius conversations, with expandable tool results and
presenter notes. Recorded mode accepts only the two saved scenarios per backend;
it does not generate new answers. Switching backends clears the conversation.
The original ECoGraph recordings remain unchanged.

### Live mode — new agent answers over cached evidence

Configure `NEBIUS_API_KEY` and `NEBIUS_MODEL` in the root `.env`, following
[.env.example](.env.example). Do not overwrite an existing working `.env`.

```dotenv
NEBIUS_API_KEY=your-nebius-api-key
NEBIUS_MODEL=zai-org/GLM-5.3-Flash
```

Select **Live agent**, choose a backend, and enter a prompt. For CRCSC, try:

> Which drugs does the model rank highest for CMS2?

Then:

> Show the highest-weight genes in its expression signature and explain the model's limitations.

Live prompts contact **Nebius Token Factory** and may incur charges. The OpenAI
Agents SDK is explicitly configured for Nebius's OpenAI-compatible **Chat
Completions** endpoint at `https://api.tokenfactory.nebius.com/v1/`.
`OPENAI_API_KEY` is not required, OpenAI tracing is disabled, and credentials stay
server-side. The LLM runs live; GNN scores remain cached.

## Choose a data backend

| | ECoGraph — original/default | CRCSC — CMS experiment |
| --- | --- | --- |
| Subtypes | GI.CIN, GI.HM-indel, GI.HM-SNV, GI.GS | CMS1, CMS2, CMS3, CMS4 |
| Molecular context | Binary mutation observations and subtype frequencies | Signed higher/lower-expression signatures |
| Cached drug model | `target_only`, seed42, frequency profile | `ppi_pathway`, seed42, signed profile |
| Reviewed run | `20261108_034016` | `20261108_173948` |
| Recorded walkthroughs | Two ECoGraph scenarios | Two CRCSC scenarios |

The two subtype systems are **not interchangeable**. CRCSC gene results contain
expression weights and discovery FDR, not mutation counts. Each agent instance
uses one fixed backend; missing data never silently falls back to the other.

Terminal commands:

```sh
# Original workflow; also the default for python -m agent:
python -m agent --backend ecograph

# CRCSC CMS workflow:
python -m agent --backend crcsc

# Check local tool data without credentials or provider calls:
python -m agent --backend ecograph --check-data
python -m agent --backend crcsc --check-data

# Validate configured provider/model startup without making an API request:
python -m agent --backend crcsc --check
```

The live evidence tools require the local PrimeKG-derived node/edge tables and
reviewed result files. Large raw/processed datasets are not included in Git.
Use `--check-data` to identify unavailable backends; recordings can still replay
saved evidence offline. See [backend details](agent/BACKENDS.md),
[agent instructions](agent/README.md), and [UI instructions](ui/README.md).

## What the agent can inspect

| Tool | Purpose |
| --- | --- |
| `get_subtypes(disease)` | List the selected backend's CRC subtype catalog |
| `get_drug_predictions(subtype, top_k)` | Read cached scores, preserving all cutoff ties |
| `get_drug_targets(drug_id)` | Retrieve direct PrimeKG targets, excluding enzyme/carrier/transporter roles |
| `get_subtype_genes(subtype, limit, gene)` | Inspect mutation observations or expression associations, depending on backend |
| `get_gene_phenotypes(gene, limit)` | Retrieve explicit recorded phenotype–gene associations |
| `get_model_evaluation()` | Retrieve the selected experiment's baseline and ablation results |

No-match results are not proof of biological absence. Graph connections and
retrieved paths are not causal explanations of model scores.

## Research findings and limitations

The CRCSC representation experiment maps **18,664 measured genes** to PrimeKG
and constructs four reproducible CMS expression signatures. Its drug-scoring
prototype adds **four subtype nodes and 400 typed expression edges** to a runtime
PrimeKG view. The full base graph is referenced read-only rather than duplicated.

The matched held-out indication results remain below popularity:

| Experiment | Molecular GNN mean AP | Drug-only ablation mean AP | Popularity AP |
| --- | ---: | ---: | ---: |
| ECoGraph CRC | 0.8254 | 0.8257 | **0.8941** |
| CRCSC CMS, PPI + pathway | 0.8250 | 0.8223 | **0.8941** |

GNN values are three-seed means; the agent's cached predictions use seed42.
**AP measures recorded-indication recovery against sampled unknown pairs—not
accuracy, efficacy probability, or subtype treatment response.** Both experiments
retain a leading group of **29 tied drugs sharing NR3C1 as their recorded target**.
The tested architecture cannot distinguish identical target neighborhoods using
those inputs alone.

Technical integration and scoring execute successfully, but useful subtype-specific
therapeutic prioritization has not been demonstrated. The CMS drug model transfers
from general disease indications without CMS therapeutic labels. Signed expression
is not drug-action direction or evidence of tumor dependency. CMS reproduction
from expression is also partly circular because CMS labels derive from expression.

A separate gene–phenotype completion benchmark found useful pathway-heuristic
signal, but its tested GNN did not outperform the strongest heuristic. The
breast-cancer notebook uses measured cell-line responses, but its saved GNN result
also trails the simple drug-by-subtype mean. Neither is validation of the CRC drug
hypotheses.

Start with these reports and guides:

- [ECoGraph–PrimeKG compatibility](reports/compatibility/compatibility_report.md)
- [Original CRC GNN findings](reports/crc_experiments/20261108_034016/crc_findings.md)
- [CRCSC molecular-profile assessment](reports/crcsc/20261108_172638/assessment.md)
- [CRCSC drug-scoring findings](reports/crcsc_drug_gnn/20261108_173948/findings.md)
- [Gene–phenotype experiment findings](reports/phenotype_gnn/20261108_155248/findings.md)
- [PrimeKG exploration and inference findings](docs/research/primekg_findings.md)

The CRCSC drug run includes an **unreviewed** literature queue. Literature review
can assess supporting, contradictory, or previously known evidence; it does not
independently verify therapeutic efficacy or establish novelty.

## Recordings and presentation

| Backend | Readable conversation | Structured trace | Review notes |
| --- | --- | --- | --- |
| ECoGraph | [Walkthrough](docs/demo/recordings/20261108_042554/walkthrough.md) | [Transcript](docs/demo/recordings/20261108_042554/transcript.json) | [Presenter notes](docs/demo/README.md) |
| CRCSC | [Walkthrough](docs/demo/recordings/20261108_181822/walkthrough.md) | [Transcript](docs/demo/recordings/20261108_181822/transcript.json) | [Prose qualifications](docs/demo/recordings/20261108_181822/REVIEW.md) |

Tool-evidence checks passed, but a recording is not a systematic agent-quality
evaluation. Actual responses are preserved; unsupported drug-class/structural
claims in the CRCSC answer are explicitly flagged as unverified.

```sh
# Offline terminal replay:
python scripts/rehearse_agent.py --replay docs/demo/recordings/20261108_181822/transcript.json

# Optional: record a NEW live CRCSC conversation (contacts Nebius):
python scripts/rehearse_agent.py --backend crcsc --live
```

See the [demo pack](docs/demo/README.md) and [pitch outline](docs/demo/pitch_outline.md).

## Repository layout

```text
agent/                         Agent runtime, provider configuration, backend adapters, tools
ui/                            Local Starlette server and buildless HTML/CSS/JavaScript
src/primekg/                   PrimeKG preprocessing package
analysis/                      Compatibility, CRCSC, and GNN experiment scripts
scripts/                       Preprocessing, initial baselines, recording, explicit live checks
notebooks/
  integration/                 Teammate PrimeKG–ECoGraph merge notebook
  models/colorectal/            Original CRC GNN notebook
  models/breast/                HER2/triple-negative response notebook
ECoGraph/                      Separate upstream project, preserved
data/                         Local datasets (see paths below)
reports/
  primekg/                     PrimeKG inventory and quality checks
  compatibility/               ECoGraph compatibility results
  inference/                   Initial gene–phenotype heuristic baselines
  crc_experiments/              Original CRC diagnostics
  crcsc/                       CMS molecular-profile reproducibility
  crcsc_drug_gnn/               CMS drug-scoring experiments and required small caches
  phenotype_gnn/               Gene–phenotype robustness benchmark
artifacts/                     Local model/prediction artifacts
docs/                         Research notes, demo recordings, repository guide
requirements/                  Agent, lightweight research, and modeling dependencies
tests/agent/                   Offline SDK, backend, UI, and recording tests
tests/research/                Analysis and model-helper tests
```

Local input locations:

- `data/raw/primekg/kg.csv`
- `data/processed/primekg/` — node and phenotype/disease/gene tables
- `data/CRCSC_data/` — uploaded expression, clinical, and CMS-label files
- `data/external/geo/GPL570.annot.gz` — official probe annotation
- `data/external/hgnc/` and `data/external/gdsc/` — modeling references/screens

Raw datasets, external references, and processed data are ignored by Git. Small
agent cache files and verified recordings are versioned. The broader research
checkout may also contain uncommitted experiment code/results; inspect `git status`
before assuming they are available in another clone or branch.

Historical manifests may retain pre-reorganization paths. See
[repository structure and migration notes](docs/repository_structure.md).

## Reproduce research

Use a separate research/modeling environment. `requirements/research.txt` contains
pandas, NumPy and SciPy; GNN/classifier workflows additionally need PyTorch,
scikit-learn and threadpoolctl. `requirements/notebooks.txt` includes the modeling
and Jupyter packages. None is required just to replay the UI's saved conversations.

Read the relevant guide before rerunning: original preprocessing/compatibility
commands can overwrite derived reports, while newer experiment runners create
new timestamped output directories.

| Workflow | Guide or entry point |
| --- | --- |
| PrimeKG preparation | `python3 scripts/prepare_primekg.py` |
| Initial annotation baselines | `python3 scripts/explore_inference.py` |
| ECoGraph compatibility | [Analysis guide](analysis/README.md) |
| Original CRC experiments | [CRC experiment guide](analysis/CRC_EXPERIMENTS.md) |
| CRCSC profiles | [CRCSC experiment guide](analysis/CRCSC_EXPERIMENT.md) |
| CRCSC drug scoring | [CMS drug GNN guide](analysis/CRCSC_DRUG_GNN.md) |
| Gene–phenotype robustness | [Phenotype GNN guide](analysis/PHENOTYPE_GNN.md) |
| Teammate notebooks | [Notebook workflow](notebooks/README.md) |

The original notebooks and historical outputs were preserved. Scripted experiments
are separate from notebook execution; no research model is deployed as a live
inference service.

## Verification

```sh
# In biohack-agent; no paid provider calls:
python -m unittest discover -s tests/agent -p 'test_*.py' -v
node --check ui/static/app.js

# In the research/modeling environment:
python3 -m unittest discover -s tests/research -p 'test_*.py' -v
```

Last verified: **26 agent/UI tests** and **25 research tests** passed in their
respective environments. Both backend data checks pass. Browser checks exercised
all four recorded scenarios, backend selection, table rendering, and desktop/mobile
layouts. Offline tests use mocked provider calls; the saved recordings were
separately executed against live Nebius.

`python -m pytest` requires an environment with both the agent and research/modeling
dependencies. Discovery is restricted to `tests/`; explicit live smoke checks are
under `scripts/smoke/` and are not run automatically.

## Preserved original agent

Development is on **`CRCSC-branch`**. The pre-integration working agent/UI is saved
at commit **`330c383`**, referenced by **`agent-ecograph-checkpoint`**. The original
ECoGraph workflow remains available without switching branches: it is still the
default backend.

Backend integration is committed as `28c5de0`; CRCSC recorded conversations are
committed as `e45cfa2`. Save outstanding work before switching branches—uncommitted
changes can follow a checkout. The checkpoint branch remains local unless pushed separately. Pushing code
does not deploy the web application; it still runs locally.
