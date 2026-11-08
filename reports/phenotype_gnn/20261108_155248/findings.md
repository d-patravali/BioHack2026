# Gene–phenotype prediction under missing evidence

## Outcome

**The tested GNN does not exceed the strongest heuristic in mean MRR. This run does not establish a learned-model advantage.** This is an executed exploratory benchmark, not a new biological discovery or an independently validated novelty claim. The CRC app, notebooks, historical results, and datasets were not modified by this experiment.

## Target and protocol

Predict explicitly recorded `phenotype_protein` associations for genes whose **entire direct phenotype annotation set is withheld**. There are 3,330 canonical positive pairs; candidates are all 15,311 phenotype nodes. Disease paths are never promoted to positive labels.

The previous gene-holdout test masks (seeds [17, 29, 43]) are reused and checked exactly by original canonical node IDs. A further 20% of annotated genes is reserved for validation, leaving approximately 60/20/20 training/validation/test genes. Earlier results used 80/20, so their numerical values are not a directly matched comparison; all heuristics here are recomputed with the same reduced supervision as the GNN. These test masks were already inspected historically: this is exploratory evaluation, not an untouched confirmatory test.

| seed | train_genes | validation_genes | test_genes | train_pairs | validation_pairs | test_pairs | test_pairs_without_training_phenotype |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 17 | 973 | 324 | 324 | 2029 | 670 | 631 | 35 |
| 29 | 973 | 324 | 324 | 1995 | 646 | 689 | 37 |
| 43 | 973 | 324 | 324 | 2003 | 648 | 679 | 45 |

All query genes remain in the molecular graph: this is **transductive annotation-cold evaluation**, not unseen-node generalization or a test on naturally unannotated genes. Three masks reuse overlapping genes and one fixed initialization seed (42); means are descriptive, not confidence intervals.

## Model and leakage controls

A small PyTorch encoder uses typed learned ID embeddings (width 24), two relation-specific mean-message layers with residual state, and a gene–phenotype dot-product decoder with phenotype bias. Training uses balanced positive and sampled-unlabeled logistic terms (five unlabeled draws per positive). Unlabeled draws exclude training positives; they are optimization contrasts, not verified biological negatives. No held-out labels enter sampling or message passing.

- **annotations_only:** training gene–phenotype edges and their reverse messages, without PPI/pathways/disease context; a learned-ID/annotation control.
- **molecular:** adds PPI and gene–pathway relations in both directions where applicable.
- **with_disease:** also adds disease–gene and positive disease–phenotype relations. Shared underlying evidence can make this channel semantically leaky even after target removal.

Exact `(type, source, id)` aliases are canonicalized before splitting. Reciprocal rows are verified and deduplicated. Held-out target edges and reverse copies are absent from the training adjacency. Checkpoints are selected using validation MRR only, evaluated every 10 epochs up to 60. No phenotype hierarchy/name features, patient data, drug labels, or CRC subtype nodes are used.

| seed | condition | best_epoch | validation_mrr |
| --- | --- | --- | --- |
| 17 | annotations_only | 50.0000 | 0.1388 |
| 17 | molecular | 60.0000 | 0.1177 |
| 17 | with_disease | 50.0000 | 0.1000 |
| 29 | annotations_only | 60.0000 | 0.1532 |
| 29 | molecular | 60.0000 | 0.1260 |
| 29 | with_disease | 60.0000 | 0.1225 |
| 43 | annotations_only | 50.0000 | 0.1301 |
| 43 | molecular | 60.0000 | 0.1162 |
| 43 | with_disease | 50.0000 | 0.1145 |

## Matched intact-evidence results

Mean over masks; per-mask and per-pair values are saved separately. MRR averages reciprocal target rank; Hits@K is the expected fraction recovered within K ranks. Other recorded positives for each query gene are filtered **only during evaluation**. Exact score ties use the existing expected-rank implementation, avoiding arbitrary ordering. Missing associations remain unknown, so these measure annotation recovery, not biological specificity.

| method | condition | mrr | gene_macro_mrr | hits1 | hits10 | hits100 |
| --- | --- | --- | --- | --- | --- | --- |
| disease_paths | heuristic | 0.0173 | 0.0208 | 0.0019 | 0.0405 | 0.2248 |
| equal_ppi_pathway | heuristic | 0.1676 | 0.1831 | 0.0965 | 0.3297 | 0.5183 |
| phenotype_popularity | heuristic | 0.1427 | 0.1587 | 0.0664 | 0.3216 | 0.8140 |
| ppi_transfer | heuristic | 0.0969 | 0.1078 | 0.0577 | 0.1889 | 0.2422 |
| relational_gnn | annotations_only | 0.1253 | 0.1403 | 0.0450 | 0.2887 | 0.8041 |
| relational_gnn | molecular | 0.1154 | 0.1466 | 0.0346 | 0.3023 | 0.8004 |
| relational_gnn | with_disease | 0.1174 | 0.1569 | 0.0431 | 0.2866 | 0.7958 |
| shared_pathway_transfer | heuristic | 0.1686 | 0.1833 | 0.1023 | 0.3168 | 0.4769 |

The strongest comparison is against pathway transfer and fixed PPI/pathway fusion, not merely a zero-score lookup. Popularity ignores the query gene. The learned annotation-only control tests whether context adds value beyond the encoder's learned IDs and label prior.

The matched pathway baseline has mean MRR **0.1686**, compared with popularity **0.1427** and the best mean GNN condition **0.1253**. The molecular GNN and the disease-enabled GNN do not establish useful incremental graph learning in this run. Context-availability strata (PPI/pathway present versus absent) are also saved in `context_strata.csv`; they are descriptive and not matched-degree causal comparisons.

## Evidence-removal experiments

Each selected checkpoint is evaluated after removing an entire relation channel, including its reverse, **at inference without retraining or checkpoint reselection**. This tests distribution-shift sensitivity, not the causal truth or independent value of a biological mechanism. Unavailable channels produce an unchanged control. `with_disease` versus `molecular` provides a separately retrained comparison with/without disease channels, although other architectural details/parameters differ.

| condition | scenario | mrr | hits10 | hits100 |
| --- | --- | --- | --- | --- |
| annotations_only | intact | 0.1253 | 0.2887 | 0.8041 |
| annotations_only | no_disease | 0.1253 | 0.2887 | 0.8041 |
| annotations_only | no_pathways | 0.1253 | 0.2887 | 0.8041 |
| annotations_only | no_ppi | 0.1253 | 0.2887 | 0.8041 |
| annotations_only | no_training_annotations | 0.0772 | 0.1857 | 0.7825 |
| molecular | intact | 0.1154 | 0.3023 | 0.8004 |
| molecular | no_disease | 0.1154 | 0.3023 | 0.8004 |
| molecular | no_pathways | 0.1146 | 0.2863 | 0.7899 |
| molecular | no_ppi | 0.1162 | 0.3019 | 0.7983 |
| molecular | no_training_annotations | 0.0905 | 0.2378 | 0.7743 |
| with_disease | intact | 0.1174 | 0.2866 | 0.7958 |
| with_disease | no_disease | 0.1046 | 0.2390 | 0.7820 |
| with_disease | no_pathways | 0.1129 | 0.2875 | 0.7942 |
| with_disease | no_ppi | 0.1175 | 0.2911 | 0.7962 |
| with_disease | no_training_annotations | 0.1013 | 0.2264 | 0.7871 |

Heuristics receive the analogous channel removals (`metrics.csv`). Removing training annotations from the GNN leaves learned parameters intact: it is not erasing every trace of training knowledge. It is a structural sensitivity test.

## Exploratory abstention

Confidence is the unfiltered top-1 minus top-2 score margin. Thresholds are chosen from **validation genes only**, at desired retained fractions 100%, 75%, 50%, and 25%, and transferred unchanged to test genes. All threshold ties are retained, so actual test coverage can differ substantially. No score margin is called a calibrated probability. We compare accepted-gene macro MRR with all-gene macro MRR, and include identical diagnostics for the heuristics.

| condition | desired_validation_coverage | test_coverage | accepted_gene_macro_mrr | all_gene_macro_mrr |
| --- | --- | --- | --- | --- |
| annotations_only | 0.2500 | 0.2562 | 0.1404 | 0.1403 |
| annotations_only | 0.5000 | 0.4887 | 0.1397 | 0.1403 |
| annotations_only | 0.7500 | 0.7407 | 0.1440 | 0.1403 |
| annotations_only | 1.0000 | 1.0000 | 0.1403 | 0.1403 |
| disease_paths | 0.2500 | 0.3045 | 0.0149 | 0.0208 |
| disease_paths | 0.5000 | 0.5103 | 0.0169 | 0.0208 |
| disease_paths | 0.7500 | 1.0000 | 0.0208 | 0.0208 |
| disease_paths | 1.0000 | 1.0000 | 0.0208 | 0.0208 |
| equal_ppi_pathway | 0.2500 | 0.2510 | 0.2442 | 0.1831 |
| equal_ppi_pathway | 0.5000 | 0.5082 | 0.2543 | 0.1831 |
| equal_ppi_pathway | 0.7500 | 1.0000 | 0.1831 | 0.1831 |
| equal_ppi_pathway | 1.0000 | 1.0000 | 0.1831 | 0.1831 |
| molecular | 0.2500 | 0.2181 | 0.1469 | 0.1466 |
| molecular | 0.5000 | 0.4877 | 0.1490 | 0.1466 |
| molecular | 0.7500 | 0.7634 | 0.1449 | 0.1466 |
| molecular | 1.0000 | 1.0000 | 0.1466 | 0.1466 |
| phenotype_popularity | 0.2500 | 1.0000 | 0.1587 | 0.1587 |
| phenotype_popularity | 0.5000 | 1.0000 | 0.1587 | 0.1587 |
| phenotype_popularity | 0.7500 | 1.0000 | 0.1587 | 0.1587 |
| phenotype_popularity | 1.0000 | 1.0000 | 0.1587 | 0.1587 |
| ppi_transfer | 0.2500 | 0.3014 | 0.1739 | 0.1078 |
| ppi_transfer | 0.5000 | 1.0000 | 0.1078 | 0.1078 |
| ppi_transfer | 0.7500 | 1.0000 | 0.1078 | 0.1078 |
| ppi_transfer | 1.0000 | 1.0000 | 0.1078 | 0.1078 |
| shared_pathway_transfer | 0.2500 | 0.2438 | 0.3681 | 0.1833 |
| shared_pathway_transfer | 0.5000 | 0.6409 | 0.2554 | 0.1833 |
| shared_pathway_transfer | 0.7500 | 1.0000 | 0.1833 | 0.1833 |
| shared_pathway_transfer | 1.0000 | 1.0000 | 0.1833 | 0.1833 |
| with_disease | 0.2500 | 0.2377 | 0.1555 | 0.1569 |
| with_disease | 0.5000 | 0.4990 | 0.1398 | 0.1569 |
| with_disease | 0.7500 | 0.7397 | 0.1515 | 0.1569 |
| with_disease | 1.0000 | 1.0000 | 0.1569 | 0.1569 |

The GNN margin policy does not show a clear, consistent selective-quality gain across the tested retention levels. As an exploratory contrast, pathway transfer at desired 25% retention accepts about 24% of test genes on average with higher accepted-gene MRR; that is selection of an easier subset, not improved full-coverage prediction, and requires fresh validation before use.

Selective MRR is annotation-recovery quality on a selected subset, not a guarantee of correct novel associations. Smaller accepted subsets can be easier; an increase does not establish overall model improvement, error control, or agent research utility. Rows with zero accepted genes have undefined quality. Constant popularity margins cannot meaningfully select among genes. A deployable abstention policy needs fresh validation and comparison at matched actual coverage.

## Limitations and interpretation

1. Only a small, study-biased set of genes/phenotypes has direct labels. Performance on these masked genes may not transfer to naturally unannotated genes; some test phenotypes have no training labels.
2. PrimeKG lacks per-edge citations/dates needed to establish independent-source or temporal validation. Disease context may duplicate the knowledge behind direct targets. PPI/pathway-only evaluation is a useful restriction, not proof of source independence.
3. Learned ID embeddings are transductive and may overfit sparse labels. One width, optimizer, and modest epoch budget do not exhaust GNN designs; poor results do not establish that GNNs are impossible.
4. Withheld links are explicit graph annotations, not experimentally verified causal effects. No prediction here establishes patient diagnosis, variant pathogenicity, drug efficacy, or novel biology.
5. Baseline votes, sparse label coverage, degree effects, and score priors can explain apparent gains. Degree-preserving graph nulls, stronger non-GNN learned fusion, additional initialization seeds, and a fresh independent test remain future controls.
6. Relation removal is a controlled whole-channel stress test, not a realistic model of every missing-data mechanism. Confidence margins and post-hoc channel sensitivity are not faithful causal explanations.
7. PhenoLinker and other work already study phenotype–gene prediction. Any claim of novelty requires a focused comparison; the present contribution is an evaluated prototype and robustness case study.

## Recommendation

Use these measured comparisons to decide whether a new GNN backend merits integration into the agent. Do not replace the current backend or present this as validated phenotype-aware inference merely because training executed. If heuristics remain stronger, keep them as the reference and focus the agent on evidence retrieval and measured uncertainty; consider stronger supervision/features only with a new validation plan. If a learned gain appears, require paired per-mask consistency and further controls before claiming generalization.

The next milestone is a **fresh, source-aware benchmark** and a separately evaluated agent investigation, not clinical deployment. This experiment deliberately leaves the working CRC demonstration unchanged.

## Reproduce and inspect

```sh
python3 analysis/phenotype_gnn.py
python3 analysis/phenotype_report.py reports/phenotype_gnn/20261108_155248
python3 -m unittest discover -s tests/research -p 'test_phenotype_gnn.py' -v
```

Dependencies: existing research requirements plus PyTorch (CPU). No PyTorch Geometric, external API, downloads, or agent credentials required. The runner creates a new timestamped folder; an existing explicit `--output` is rejected. The report generator writes only this run's `findings.md` and summary tables.

Artifacts: `metrics.csv`, `per_pair_metrics.csv.gz`, `splits.csv`, `training_curves.csv`, `selective_prediction.csv`, `top_hypotheses.csv`, `context_strata.csv`, `model_*.pt`, and `manifest.json`. Hypotheses are bounded top-five examples; arbitrary within-tie order is not meaningful. The manifest records settings, versions, input/code hashes and preservation checks. Runtime: 139.9 seconds. Protected input hashes unchanged: **True**.
