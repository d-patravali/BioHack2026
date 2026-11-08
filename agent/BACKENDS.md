# Selectable research backends

The original ECoGraph agent remains the default. CRCSC is an explicit alternative;
changing the backend does not retrain a model or replace the original data.

## Run

```sh
conda activate biohack-agent
# Original GI genomic-subtype workflow:
python -m agent --backend ecograph --check-data
python -m agent --backend ecograph
# New expression-defined CMS workflow:
python -m agent --backend crcsc --check-data
python -m agent --backend crcsc
# Local UI with both choices:
python -m ui
```

Open http://127.0.0.1:8765. Choose a **Data backend** in either mode.
Backend changes clear conversation history and update the subtype catalog,
guided prompts, and evaluation panel. Recorded mode replays the saved conversation for the selected backend without
network access. CRCSC has its own actual Nebius recording, separate from the
original ECoGraph recording; neither is presented as live inference.

Sample CRCSC prompt:

> Which drugs does the model rank highest for CMS2?

Follow-up:

> Show the highest-weight genes in its expression signature and explain the model's limitations.

| Property | ECoGraph (default) | CRCSC |
| --- | --- | --- |
| Catalog | GI.CIN, GI.HM-indel, GI.HM-SNV, GI.GS | CMS1, CMS2, CMS3, CMS4 |
| Molecular evidence | Mutation observation frequencies and patient counts | Higher/lower relative expression, signed weights, discovery FDR |
| Reviewed run | `reports/crc_experiments/20261108_034016` | `reports/crcsc_drug_gnn/20261108_173948` |
| Cached model | seed42, target_only, frequency | seed42, ppi_pathway, signed |
| Saved offline walkthroughs | Two original scenarios | Two CRCSC scenarios |

All six tool names are retained. Each agent instance binds tools to a fixed
backend; there is no mutable process-wide selector. Wrong-catalog subtype names
return unsupported_subtype; missing CRCSC files return backend_unavailable,
never silently falling back to ECoGraph. Shared targets and phenotype annotations
come from the same read-only PrimeKG tables. An HTTP conversation cannot change
backend without a reset.

## Preserved working version

Before integration, the existing UI and its passing tests were checkpointed:

- Commit **`330c383`** — `Checkpoint working ECoGraph agent UI before CRCSC integration`.
- Local branch **`agent-ecograph-checkpoint`** points to that commit.
- Current development remains on **`CRCSC-branch`**.

No checkout is necessary to use the old workflow: leave the default selection or
pass `--backend ecograph`. The Git reference also preserves the pre-integration
code; save other work before any branch checkout. Raw datasets remain local
ignored inputs and are not packaged into a commit.

## Scientific boundaries

Both backends serve **cached experimental scores**, preserving full cutoff ties
(29 leading compounds in the reviewed snapshots). Raw logits are not response
probabilities. Neither drug GNN used phenotype edges. No literature review or
therapeutic validation is implied.

CRCSC genes are expression associations, not mutation frequencies or patient-gene
mutation counts. Signed pooling transfers a decoder trained on unsigned general
disease associations to expression-derived subtypes. That transfer, action
direction, and therapeutic benefit remain unvalidated. The current model does
not establish useful subtype-specific leading candidates; popularity still wins
on its sampled held-out indication benchmark. The UI shows the selected run's
three-seed mean; its cached predictions use seed42, whose AP differs from that mean.

## Files and verification

Implementation: `agent/backends.py`, backend-bound tool factory in `agent/tools.py`,
`create_agent(..., backend=...)`, CLI `--backend`, and the local UI selector.
The CRCSC adapter needs only standard-library CSV/gzip/JSON; PyTorch is not needed
to use the agent. Minimal reviewed cache files are versioned with this integration;
large model checkpoints and source expression data are not required by the agent.
The usual PrimeKG node/evidence tables are still required locally, as before.

```sh
python -m unittest discover -s tests/agent -p 'test_*.py' -v
node --check ui/static/app.js
```

Tests cover both catalogs, fixed source/model selection, ties, expression semantics,
missing-data behavior, per-agent tool isolation, Nebius Chat Completions routing,
HTTP history isolation, and original recorded replay. Provider calls are mocked;
these checks do not establish agent scientific truthfulness or drug efficacy.

## CRCSC recorded conversation

Select **Recorded**, choose **CRCSC · CMS expression**, then click either guided
investigation. The saved session is at
`docs/demo/recordings/20261108_181822/transcript.json`; read its `REVIEW.md` for
prose qualifications. Tool evidence checks passed; all graph/results hashes stayed
unchanged. Drug-class/structural-similarity descriptions in the saved answer were
not supplied by tools and remain unverified; the replay flags them explicitly.

To record a new live session (uses Nebius and may incur charges):

```sh
python scripts/rehearse_agent.py --backend crcsc --live
```

Recordings are pinned by backend in the UI. An incomplete or wrong-backend
recording is unavailable, never silently replaced by the other dataset's session.
