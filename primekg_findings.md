# Does PrimeKG support a GNN-centered biological inference project?

**Recommendation:** pursue a small, falsifiable **gene–phenotype annotation-completion benchmark**, with a pathway-aware heterogeneous GNN as a possible next experiment. PrimeKG contains measurable predictive structure, but its direct phenotype labels are sparse and biased. The current evidence supports a benchmark prototype, not a claim of new biological discovery or GNN superiority.

## Context and scope

The hackathon asks: **“How can we prototype agentic frameworks that leverage GNNs to parse phenotype-genotype knowledge graphs?”** Our working direction is biological hypothesis investigation rather than patient diagnosis. We want a learned graph model to earn its place by improving prediction over strong simple methods; an agent could subsequently orchestrate queries and evidence retrieval.

[SHEPHERD](https://pmc.ncbi.nlm.nih.gov/articles/PMC12181314/) already studies rare-disease diagnosis using an enriched knowledge graph, simulated patients, and geometric deep learning. It is relevant prior work, but this experiment neither reproduces nor compares against its patient-level tasks. Our question is narrower: **can indirect biological connections recover explicitly recorded gene–phenotype links when those links are hidden?**

All numerical findings below come from the supplied `kg.csv` and executed experiments. External sources are used only to interpret construction/provenance and prior work. No GNN or agent was trained or implemented.

## 1. What is actually in this CSV?

There are **8,100,498 ordered edge rows**, **129,375 node indices**, **10 node types**, and **30 distinct `relation` labels**. Including endpoint types and `display_relation` produces 56 relation/display/type combinations. Counts were checked by rescanning the raw CSV against the existing preprocessing inventory.

| Node type | Nodes |
| --- | ---: |
| `biological_process` | 28,642 |
| `gene/protein` | 27,671 |
| `disease` | 17,080 |
| `effect/phenotype` | 15,311 |
| `anatomy` | 14,035 |
| `molecular_function` | 11,169 |
| `drug` | 7,957 |
| `cellular_component` | 4,176 |
| `pathway` | 2,516 |
| `exposure` | 818 |

Each row has `relation`, `display_relation`, and five fields per endpoint: `{x,y}_index`, `{x,y}_id`, `{x,y}_type`, `{x,y}_name`, `{x,y}_source`. Endpoint indices identify nodes within this file. External IDs are opaque strings, scoped by type and source; names are display metadata. `source` identifies the node namespace, not the publication supporting an edge. The CSV has no per-edge dates, citations, or confidence values.

Gene nodes use `NCBI`, phenotypes use `HPO`, and diseases use `MONDO` (15,813 nodes) or `MONDO_grouped` (1,267). Genes and proteins share one node type. The phenotype type also includes drug side effects, so it should not automatically be interpreted as a collection of rare-disease symptoms.

The existing pipeline found no blank fields, raw self-loops, or conflicting metadata for an index. However, **122 gene nodes represent 61 repeated external IDs**, and four anatomy nodes also share IDs. For experiments only, exact `(type, source, id)` matches are mapped to their lowest original index. This yields **27,610 canonical genes**, 15,311 phenotypes, and 17,080 diseases. Original preprocessing outputs and raw data remain unchanged. No name-based or cross-namespace equivalences are assumed.

### Relevant edge families

The following counts are distinct canonical pairs. Reciprocal CSV rows count as one association for this table. Same-type pairs are unordered here only to describe storage; this does not establish biological direction.

| Relation (`display_relation`) | Distinct pairs | Why it matters |
| --- | ---: | --- |
| `phenotype_protein` (`associated with`) | **3,330** | Explicitly recorded target labels |
| `disease_protein` (`associated with`) | 80,411 | Gene–disease context |
| `disease_phenotype_positive` (`phenotype present`) | 150,317 | Positive disease–phenotype context |
| `disease_phenotype_negative` (`phenotype absent`) | 1,193 | Distinct negative annotations; not patient observations |
| `protein_protein` (`ppi`) | 321,075 | Transfer of phenotype evidence through interacting genes |
| `pathway_protein` (`interacts with`) | 42,646 | Transfer through shared pathways |

The six families above have verified reciprocal rows. For example, 3,330 direct pairs occupy 6,660 raw rows. Across the earlier three-type extraction, all 1,214,512 ordered edges have reciprocal copies and none is an exact ordered duplicate.

Other available relations, with raw counts counting both endpoint orders where present:

| Family | Exact relation labels and raw rows |
| --- | --- |
| Functional annotations | `bioprocess_protein` 289,610; `molfunc_protein` 139,060; `cellcomp_protein` 166,804 |
| Hierarchies (`parent-child` display label) | `disease_disease` 64,388; `phenotype_phenotype` 37,472; `pathway_pathway` 5,070; `bioprocess_bioprocess` 105,772; `molfunc_molfunc` 27,148; `cellcomp_cellcomp` 9,690; `anatomy_anatomy` 28,064; `exposure_exposure` 4,140 |
| Anatomy expression | `anatomy_protein_present` 3,036,406; `anatomy_protein_absent` 39,774 |
| Drug context | `drug_drug` 2,672,628; `drug_protein` 51,306; `drug_effect` 129,568; `indication` 18,776; `contraindication` 61,350; `off-label use` 5,136 |
| Exposure context | `exposure_bioprocess` 3,250; `exposure_cellcomp` 20; `exposure_disease` 4,608; `exposure_molfunc` 90; `exposure_protein` 2,424 |

`drug_protein` has distinct display labels (`target`, `enzyme`, `carrier`, `transporter`); retaining only `relation` would discard those distinctions. Full typed counts are in `reports/relationships.csv`. PPI, pathways, GO functions/processes, and tissue expression offer molecular context beyond phenotype matching. Whether each helps is an experimental question; only PPI and pathway transfer were tested here.

## 2. Direct annotations are present, but very sparse

“Direct” means an explicit `phenotype_protein` CSV edge—not a causal assertion, an experimental validation, or a claim that the underlying publication measured that relationship directly. The upstream [PrimeKG construction notebook](https://github.com/mims-harvard/PrimeKG/blob/main/knowledge_graph/build_graph.ipynb) builds these links from DisGeNET records categorized as phenotypes, mapped to HPO. This gives construction context, but the local CSV cannot identify the supporting record for each edge.

| Coverage measure | Observed value |
| --- | ---: |
| Genes with direct phenotype annotations | **1,621 / 27,610 (5.87%)** |
| Phenotypes appearing in direct annotations | **356 / 15,311 (2.33%)** |
| Direct annotations per annotated gene | Median **1**, maximum **27** |
| Annotated genes with exactly one direct phenotype | **958** |
| Genes with gene–disease–positive-phenotype paths | **8,532 (30.90%)** |
| Distinct gene–phenotype pairs reached by these paths | **820,900** |
| Genes reached by paths but without direct annotations | **7,086** |
| Genes with either direct or indirect coverage | **8,707 (31.54%)** |
| Directly annotated genes without these disease paths | **175** |
| Direct pairs also reachable through disease paths | **1,165 / 3,330 (34.98%)** |

Most genes have neither kind of coverage. The indirect expansion is much larger than the target set, but **819,735 path-derived pairs lack a recorded direct target edge**. These are unknown associations, not 819,735 new discoveries and not validated training labels.

Molecular context has wider coverage: 18,354 canonical genes have PPI neighbors, 10,849 have pathway annotations, and 9,080 have disease associations. This motivates testing transfer from molecular context for genes whose direct annotations are masked.

Concrete examples, using names exactly as stored:

- **Explicit edge plus path:** `SORT1` (index 9) has a direct link to `Elevated hepatic transaminase` (24475). It also connects through `benign recurrent intrahepatic cholestasis` (31033) to that phenotype.
- **Path without direct edge:** `ZRSR2` (2) connects through `myelomonocytic leukemia` (28502) to `Abnormal bleeding` (22180), but has no direct `phenotype_protein` edge to that phenotype in this file. This is a hypothesis-generating connection, not a documented gene–phenotype label.

Some direct labels are common: `Seizure` has 218 associated genes, `Hearing impairment` 115, and `Spontaneous abortion` 109. The limited phenotype support and uneven frequencies make a phenotype-popularity baseline essential.

## 3. Baseline experiment

### Target and evaluation design

Given a gene, rank **all 15,311 phenotype nodes** to recover hidden, explicitly recorded `phenotype_protein` pairs. No indirect path is used as a positive label. No missing edge is declared biologically negative.

Two reproducible 80/20-style masking schemes were evaluated with seeds **17, 29, 43**:

1. **Edge holdout:** hide 666 of the 3,330 unique canonical pairs; train on 2,664. There are 508–525 evaluated genes per seed. This tests annotation completion, allowing some query genes to retain other labels.
2. **Gene annotation holdout:** select 324 of the 1,621 annotated genes and hide **all** their direct phenotype pairs. This produces 631, 689, and 679 test pairs. Their PPI, pathway, and disease context stays available. This is annotation-cold evaluation, not generalization to unseen graph nodes or a representative sample of naturally unannotated genes.

The graph is represented by canonical unique pair sets before splitting, eliminating reciprocal and exact-ID alias leakage. All direct-label votes and phenotype frequencies use training pairs only. Held-out labels cannot reappear through their reverse row or an exact-ID duplicate. Shared-pathway scores exclude the query gene's own votes; PPI scores exclude itself. Hierarchy edges and phenotype-name features are not used, avoiding ontology propagation shortcuts in these baselines.

Gene–disease–phenotype paths remain intentionally available for the disease-path method. They are the signal under test, but may encode overlapping source knowledge; independence of the underlying evidence cannot be established. PPI/pathway methods exclude disease paths entirely, providing an important ablation.

### Methods

| Method | Score for candidate phenotype p and query gene g |
| --- | --- |
| Direct lookup | 1 for a training direct annotation, otherwise 0; a retrieval-only control |
| Phenotype popularity | Number of training genes directly annotated with p; ignores g |
| Disease paths | Sum over diseases linked to g and p of `1 / number of positive phenotypes of that disease` |
| PPI transfer | Number of PPI neighbors of g with a training direct annotation to p |
| Shared-pathway transfer | For each pathway of g, count other member genes annotated with p, divided by the number of other members; sum across pathways |
| Equal PPI/pathway mix | Normalize each query's PPI and pathway vectors by its channel maximum, then average with fixed equal weights |

There is no parameter tuning or learned model. All methods score the same candidate universe. Evaluation is **filtered**: other known direct positives for the query gene are excluded from ranking, retaining the test target. Held-out positives are used for this evaluation filter only, not for prediction. Unannotated candidates may include true associations, so this measures recovery of existing annotations, not biological specificity or precision.

**MRR** averages the reciprocal rank of each held-out pair; higher means known annotations are placed nearer the top. **Hits@K** is the fraction of targets in the top K. Tied scores use the exact expectation under uniformly random tie ordering, preventing arbitrary node indices from deciding success. The tiny nonzero direct-lookup scores come entirely from random zero-score ties: lookup provides **no positive score for any held-out target** and would return no verified answer if it abstained.

### Results

Values are means of the three seed-level evaluations, weighted equally across seeds. Within each seed the main metrics average over held-out pairs.

| Method | Edge MRR | Edge Hits@10 | Edge Hits@100 | Gene-holdout MRR | Gene-holdout Hits@10 | Gene-holdout Hits@100 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Direct lookup / zero-score ties | 0.0007 | 0.07% | 0.65% | 0.0007 | 0.07% | 0.65% |
| Phenotype popularity | 0.1438 | 32.38% | **81.89%** | 0.1416 | 32.30% | **82.28%** |
| Disease paths | 0.0203 | 4.45% | 23.65% | 0.0173 | 4.05% | 22.48% |
| PPI transfer | 0.1127 | 22.25% | 28.81% | 0.1125 | 21.22% | 27.89% |
| Shared-pathway transfer | **0.1821** | 31.75% | 51.66% | **0.1815** | 33.00% | 51.47% |
| Equal PPI/pathway mix | 0.1813 | **33.62%** | 56.08% | 0.1800 | **35.02%** | 56.12% |

Shared-pathway MRR exceeds popularity in each paired seed. Its MRR ranges are 0.1650–0.1997 for edge holdout and 0.1770–0.1866 for gene holdout. Gene-macro MRR (equal weight per evaluated gene within each seed) is also higher for pathways than popularity: **0.1910 vs 0.1547** in edge holdout and **0.1970 vs 0.1610** in gene holdout. These repeated masks overlap; their ranges are not confidence intervals or a significance test.

Interpretation:

- **Useful structure exists:** pathway transfer improves early ranking, including when all direct annotations of query genes are removed. This is a concrete reason to investigate pathway-aware learning.
- **Coverage remains a bottleneck:** in gene holdout, nonzero scores reach only 51.21% of target pairs for pathways, 27.46% for PPI, and 34.50% for disease paths. Popularity reaches 95.01%. A graph method can be specific near the top while missing many labels entirely.
- **Broad recall favors popularity:** its Hits@100 is substantially better. A blanket claim that graph methods outperform simple baselines would be false.
- **Disease traversal alone is insufficient:** its large reachable set produces weak ranks. Reachability should not be confused with useful discrimination.
- **More relations do not automatically help:** the equal PPI/pathway mix has slightly lower MRR than pathways alone, though better Hits@10/100. A GNN would need to learn useful channel weighting rather than merely aggregate everything.

## 4. Limitations that could change the conclusion

1. **Labels are sparse and selected by prior research.** Only 356 phenotypes and 1,621 genes have direct labels; 958 genes have just one. Good results may reflect annotation practices and well-studied pathways. Performance on naturally unannotated genes cannot be measured using this file alone. Between 25 and 34 test phenotype identities per split have no training direct labels, constraining transfer methods.
2. **Evidence independence is unresolved.** Disease and direct phenotype links may share source studies or clinical knowledge. Removing explicit target edges does not eliminate all semantic leakage. Exact-ID aliases were handled, but cross-ontology synonyms, phenotype near-equivalents, and disease grouping can still encode overlap. Hierarchy-based models would require additional family/ontology-aware splits.
3. **Positive/negative disease annotations conflict.** **172 disease–phenotype pairs occur in both positive and negative relations.** Neither relation is silently discarded or reconciled. The exploratory disease-path score uses positive edges even for these conflicts and ignores negatives; it is not a clinically coherent reasoning model. Source-level review is needed before using these annotations as supervision.
4. **Missing is unknown.** Ranking distractors are unannotated, not experimentally disproven associations. MRR/Hits measure annotation recovery; they cannot validate novel biology. Exact-HPO ranking also penalizes biologically related terms without assessing semantic similarity.
5. **The split is transductive and retrospective.** Query genes remain in the context graph. No temporal, external-source, unseen-node, or disease-family generalization claim is supported. There is no variant information, inheritance model, or patient cohort here.
6. **No GNN comparison yet.** The experiment establishes a useful heuristic baseline, not a need for deep learning. Popularity, pathway weighting, and simple score fusion may capture much of the available signal. Degree-matched graph randomization, phenotype-macro metrics, and stratification by annotation frequency remain useful controls.
7. **Version/provenance uncertainty.** The local schema does not include the additional OMIM relation families described in later [upstream updates](https://github.com/mims-harvard/PrimeKG). The exact local SHA-256 is recorded in `coverage_summary.json`; current upstream documentation is not proof of this file's release date. The [PrimeKG paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9893183/) describes source integration, but the edge CSV omits the lineage needed for strict source-disjoint evaluation.

## 5. Where a GNN could reasonably add value

| Candidate task | Rationale | Required evidence before claiming success |
| --- | --- | --- |
| **Gene–phenotype annotation completion, including masked-gene annotations** | Most aligned with our question; pathway baselines already show early-ranking signal | Beat popularity, normalized pathway transfer, and simple learned score fusion on a separate validation/test design; show coverage and frequency strata |
| Disease–gene link prediction | 80,411 documented pairs offer more supervision; PPI, pathways, GO, and tissue context could explain associations | Separate benchmark with target/reverse removal and disease-family/source controls; cannot borrow the current results as evidence of performance |
| Disease–phenotype completion | 150,317 positive pairs offer a larger target set | Resolve 172 positive/negative conflicts and control hierarchy/semantic shortcuts; less directly a genotype-inference contribution |

A heterogeneous GNN might combine molecular paths, distinguish relation usefulness, and propagate information beyond fixed path templates. That is a hypothesis. Attention weights or retrieved paths would not by themselves establish causal mechanisms or faithful explanations.

## 6. Recommended next steps

1. **Adopt a conditional benchmark goal:** “Does relation-aware learning improve early gene–phenotype annotation recovery beyond pathway transfer and popularity?” Use gene annotation holdout as the primary stress test, and report both early precision-oriented ranking and broad recall. Retain the modest target scope: this is not a rare-disease diagnosis system.
2. **Audit labels before enlarging the model.** Investigate provenance, repeated IDs, disease/phenotype overlaps, annotation frequency, and positive/negative conflicts. Decide whether an external or newer independently sourced gene–phenotype label set can support stronger validation. Do not create positives from disease paths.
3. **Freeze a new validation/test protocol before model tuning.** These masks have already been inspected. Compare pathway/PPI ablations, a popularity-plus-pathway heuristic, and a small learned non-GNN fusion baseline. Add degree-preserving/null graph controls and phenotype-frequency strata. A GNN should earn a gain beyond those baselines, not just beyond direct lookup.
4. **Only then try a small relation-aware GNN.** Start with genes, phenotypes, and pathways/PPI; add GO or anatomy in measured ablations. Treat disease paths as a separately reported channel because of evidence overlap. If a GNN cannot improve a chosen metric without unacceptable coverage loss, present that result honestly or pivot to the richer disease–gene benchmark.
5. **Add an investigative agent after the benchmark is credible.** It could retrieve typed paths, run the model, compare heuristic scores, and mark evidence as explicit, indirect, or predicted. Novel hypotheses would still require independent evidence. Its role is orchestration, not the justification for the GNN.

## Reproducing and inspecting this work

Run from the repository root with the existing requirements installed:

```sh
python3 scripts/prepare_primekg.py   # only needed if processed tables are absent/stale
python3 scripts/explore_inference.py
python3 -m unittest discover -s tests
```

The analysis reuses `data/processed/nodes.csv` and the existing core extraction, rescans `kg.csv` for relation verification and pathway context, and checks target consistency. It reads raw data without modifying it. It uses pandas/NumPy, no graph framework, negative sampling, dense all-pairs graph, GNN, or LLM. Context and candidate vectors fit in memory; the raw CSV is read in chunks.

Small outputs are in `reports/inference/`: `coverage_summary.json`, `gene_coverage.csv`, `examples.json`, `baseline_metrics.csv`, `baseline_path_strata.csv`, `test_splits.csv`, and `split_audit.csv`. The manifest records test pairs; training pairs are the remaining canonical direct pairs. Versions used: pandas 2.2.2 and NumPy 2.4.4. Four small ranking tests cover filtering, ties, and unique ranks; runtime assertions verify reciprocal pairs, aliases, raw/processed consistency, and train/test separation.
