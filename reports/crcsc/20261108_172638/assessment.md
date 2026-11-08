# CRCSC assessment: a stronger molecular starting point, not validated drug inference

**Branch:** `CRCSC-branch`. **Completed run:** `20261108_172638`.

## What was established

- All 566 expression samples join to CMS labels and clinical records. The final labels provide 519 classified cases (CMS1:91, CMS2:232, CMS3:69, CMS4:127); 47 unclassified cases are excluded from supervised analysis.
- The supplied split yields 409 classified discovery and 110 classified validation cases. Features, scaling, signature selection, and classifier fitting use discovery only.
- 42,904 probes have unambiguous single-NCBI-ID mappings, yielding 20,848 genes. **18,664 map to PrimeKG (89.52%)**. This provides stable-ID bridges without fuzzy gene-name matching.
- Each CMS receives 100 discovery-selected genes with signed adjusted expression weights. Direction agreement in validation is 100–100%; **90–100 of 100** also satisfy validation BH q<=.05 in the same direction. Batch/age/sex/stage/location adjustment is descriptive and does not establish causality.
- 1,411 PrimeKG Reactome pathways have at least five measured genes. Across all pathway effects, discovery–validation Spearman correlations range **0.623–0.897**. Selected pathway direction agreement is 95–100%.

## Simple classifiers and graph controls

| Diagnostic | Validation accuracy | Balanced accuracy |
| --- | ---: | ---: |
| Majority | 42.7% | 25.0% |
| Batch only | 36.4% | 27.5% |
| Expression logistic | 91.8% | 91.7% |
| Signed 100-gene signature scores | 89.1% | 92.5% |
| Pathway logistic | 94.5% | 95.2% |

Shuffled-pathway accuracy ranges **90.0–93.6%**. The pathway-versus-expression accuracy difference is 2.73 percentage points; its paired within-class bootstrap 95% interval is **[-0.91, 6.36] percentage points**, including zero. This panel does not establish a reliable advantage from curated pathway membership. The bootstrap captures uncertainty within this validation set, not upstream label uncertainty or external-cohort generalization.

CMS labels were derived from expression, and GSE39582 contributed to the original CMS work. High reproduction accuracy is expected to some degree and must not be pitched as new subtype discovery, independent biological validation, or a demonstrated need for a GNN.

## Connection to the intended subtype–drug goal

The goal remains a model that uses subtype molecular context in a larger graph to score **subtype–drug hypotheses**. This experiment tests the reliability and compatibility of that context before reusing it in the model. It does not train a GNN or drug-response predictor.

The 100-gene signatures have the following recorded direct-target connectivity:

| Subtype | Signature genes with direct drug targets | Unique connected drugs |
| --- | ---: | ---: |
| CMS1 | 21/100 | 74 |
| CMS2 | 18/100 | 117 |
| CMS3 | 21/100 | 65 |
| CMS4 | 10/100 | 81 |

These counts are not rankings or therapeutic evidence. A drug targeting an over- or underexpressed gene may have an unsuitable mechanism, direction, tissue exposure, or toxicity profile. Expression does not establish tumor dependency. Negative expression-profile weights must not be silently turned into unsigned disease-association strengths.

## Verdict and next steps

**Use CRCSC as the preferred new CMS molecular-profile benchmark.** It is more appropriate than broad mutation unions for investigating expression-defined CMS biology, but is a complementary modality/classification, not proof that ECoGraph's genomic labels are invalid. Do not equate CMS1–4 with the four GI genomic subtypes.

The existing GNN's popularity-dominated supervision and identical drug-target neighborhoods are still unresolved. Replacing input profiles alone cannot demonstrate therapeutic prediction. Before a model extension:

1. Freeze this representation and a new evaluation plan rather than tuning repeatedly on these validation cases.
2. Check a separate expression cohort, upstream annotation/preprocessing, tissue composition, and stability under alternative gene signatures.
3. Seek independent subtype-sensitive drug-response or dependency evidence with defensible CMS labels; evaluate known/unknown associations carefully.
4. Add subtype-to-gene expression-context edges without conflating them with causal disease associations; preserve sign, subtype definition, source, and membership confidence.
5. If a GNN is then trained, compare against matched expression/pathway and drug baselines, include relation/profile ablations, and distinguish compound-specific information from shared targets.
6. Use literature review to audit supporting, contradictory, and previously known evidence with citations. Literature review can support or weaken a hypothesis, but does not independently verify efficacy; lack of a graph edge or a paper is not proof of novelty.

The new experiment did not modify the agent or original data/results, and no clinical efficacy is claimed. See [full methods/results](findings.md), [input provenance](manifest.json), and the repository guide `analysis/CRCSC_EXPERIMENT.md`.
