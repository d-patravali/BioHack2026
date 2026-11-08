# First-hour milestone: connect the research tools

The agent can now inspect cached model outputs and real recorded graph evidence.
It is not a live GNN service, and its phenotype tool does not change the trained
model. All dataset, notebook and experiment files remain read-only.

## Start and verify

From the project root:

```sh
conda activate biohack-agent
python -m agent --check-data
python -m agent --check
python -m unittest discover -s tests/agent -p 'test_*.py' -v
python -m agent
```

The data check needs no API credentials or network. Prompts contact Nebius.
The existing `NEBIUS_API_KEY` and `NEBIUS_MODEL` are read from the root `.env`.
Tool calls and returned evidence are displayed during execution.

## First walkthrough

Ask:

> Investigate the cached leading candidates for GI.CIN. Retrieve the top five,
> inspect shared drug targets, check the target gene's phenotype evidence and
> subtype mutation observations, and compare the model with popularity and the
> molecular-context ablation. What can we actually conclude?

Expected factual checks, not a script the model should invent:

- The cached seed42, target-only, frequency-profile results contain a 29-drug
  leading tie. No compound is uniquely best; requesting five includes all ties.
- Querying a member such as DB00596 (Ulobetasol) returns its recorded target
  NR3C1 / NCBI2908. It does not establish a therapeutic action or benefit.
- In this local snapshot NR3C1 has no direct `phenotype_protein` association.
  The GI.CIN retained positive-mutation profile also has no match for it.
  Neither observation proves biological absence.
- The selected seed42 AP is 0.8158, the target-only three-seed mean is 0.8254,
  drug-only mean is 0.8257, and popularity AP is 0.8941. These are sampled
  disease-indication recovery metrics, not subtype response validation.
- The model did not use phenotype edges. Retrieved graph paths/associations
  must not be claimed as model attributions or independent clinical evidence.

For a second, positive phenotype–genotype evidence example, ask:

> Retrieve the three most frequent observed genes for GI.CIN, then look up
> recorded phenotype associations for APC (NCBI324). Separate cohort observations
> from graph associations and explain whether the CRC GNN used those phenotypes.

The local snapshot has four direct APC phenotype associations, including
Gastrointestinal hemorrhage and Hematochezia. This is graph retrieval, not a new
prediction by the CRC GNN or an independent explanation of drug efficacy.

## Verification

The first walkthrough was executed against live Nebius GLM-5.3-Flash after
implementation. The agent invoked all six tools, retrieved the 29-drug tie,
checked NR3C1 evidence, and reported the popularity/ablation limitations without
claiming efficacy. Thirteen offline agent tests passed, including a complete
mocked multi-tool SDK loop. This is one successful smoke test, not a systematic
agent-quality evaluation.

Both walkthroughs now have a [recorded backup and presenter review notes](../docs/demo/README.md).
Use [the pitch outline](../docs/demo/pitch_outline.md) to prepare the presentation.

## What the tools read

| Tool | Backend |
| --- | --- |
| `get_subtypes(disease)` | Static four-subtype catalog |
| `get_drug_predictions(subtype, top_k)` | Reviewed tie-aware cached experiment CSV; cutoff 1–20, may return extra tied candidates |
| `get_drug_targets(drug_id)` | PrimeKG-derived direct-target table, filtering display role `target` |
| `get_subtype_genes(subtype, limit, gene)` | Frequency-profile observations with sample numerator/denominator; optional exact gene filter |
| `get_gene_phenotypes(gene, limit)` | Explicit PrimeKG phenotype–gene edges and node metadata; reciprocal copies deduplicated |
| `get_model_evaluation()` | Actual fixed-run AP and ablation metrics |

Later services can replace these backend adapters while retaining the tool
contracts. The next milestone is a reliable pitch walkthrough and evidence
display—not more GNN training or a claim that this tool layer improves predictions.
