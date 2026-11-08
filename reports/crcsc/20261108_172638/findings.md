# CRCSC CMS molecular-profile experiment

## Scope and outcome

This new experiment uses the uploaded GSE39582 expression/clinical/CMS files. It evaluates whether **discovery-derived gene and pathway profiles reproduce in the supplied validation group** and whether known CMS assignments can be reproduced with simple models. It does not train a GNN, merge graphs, rank drugs, or validate drug–subtype links. The existing CRC agent remains unchanged.

## Cohort and identifiers

All 566 expression samples join one-to-one to clinical/CMS records by GSM ID. Expression has 54,675 unique Affymetrix probes with finite numeric values. The final combined CMS column is used; NOLBL samples are excluded from supervised profiles/classification, not forced into a class. Network-only and RF-only labels remain in the sample audit.

| CMS | discovery | validation |
| --- | --- | --- |
| CMS1 | 76 | 15 |
| CMS2 | 185 | 47 |
| CMS3 | 52 | 17 |
| CMS4 | 96 | 31 |
| NOLBL | 34 | 13 |

Probe annotation: [official GPL570 annotation](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz); metadata dated Aug 09 2016. Accept only a single distinct numeric NCBI ID, excluding ambiguous multi-gene probes and unmapped/control probes. Aggregate multiple accepted probes per gene by within-sample median. The snapshot is old and mappings are provisional; no alias/name guesses or silent remapping of retired IDs are used. Of 20,848 annotated unique genes, **18,664 map by exact NCBI ID to PrimeKG (89.52%)**; 18,664 have nonzero discovery variance. Exact-ID PrimeKG node aliases are canonicalized. The expression header omits the probe column label; the reader explicitly handles the extra first field.

## Discovery-only profile construction

Discovery gene means/standard deviations are used to standardize both groups. For each CMS versus the other three, OLS estimates the expression contrast controlling batch, age, gender, stage, and tumor location, separately within each cohort. BH correction is per subtype and feature level. Up to 100 discovery q<=.05 genes, ranked by absolute adjusted effect, form a signed, L1-normalized profile. No validation outcome selects features, weights, or thresholds. Standardization reflects relative cohort expression, **not tumor-versus-normal expression**.

PrimeKG `pathway_protein` pairs define Reactome membership. Pathways with >=5 measured variable genes are retained (1,411). Features are mean discovery-standardized expression of measured members, not validated activation scores. Overlapping and broad pathways are correlated. Gene permutations preserve pathway sizes and membership structure but do not fix each named gene's degree.

## Profile replication

Spearman compares adjusted effects across all measured features. Sign agreement and validation BH support refer to discovery-selected signatures only; validation labels are used for evaluation, not training. Shared expression-based labeling and processing can make this easier than genuinely independent validation.

| subtype | level | tested_features | discovery_q05_features | selected_features | all_effect_spearman | selected_sign_agreement | selected_validation_q05_same_direction | discovery_residual_df | validation_residual_df |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CMS1 | gene | 18664 | 8365 | 100 | 0.7612 | 1.0000 | 90 | 392 | 96 |
| CMS1 | pathway | 1411 | 962 | 100 | 0.8966 | 1.0000 | 86 | 392 | 96 |
| CMS2 | gene | 18664 | 8662 | 100 | 0.6011 | 1.0000 | 97 | 392 | 96 |
| CMS2 | pathway | 1411 | 819 | 100 | 0.7715 | 0.9700 | 88 | 392 | 96 |
| CMS3 | gene | 18664 | 5616 | 100 | 0.6854 | 1.0000 | 100 | 392 | 96 |
| CMS3 | pathway | 1411 | 553 | 100 | 0.6232 | 0.9500 | 59 | 392 | 96 |
| CMS4 | gene | 18664 | 9486 | 100 | 0.6453 | 1.0000 | 100 | 392 | 96 |
| CMS4 | pathway | 1411 | 1013 | 100 | 0.8830 | 1.0000 | 98 | 392 | 96 |

## Prespecified subtype-reproduction diagnostics

Fit on labeled discovery samples, evaluate on labeled validation samples. No test-tuned hyperparameters. Gene logistic/nearest-centroid models use the top 2,000 discovery-variance genes. Logistic models have C=1, class-balanced loss, and max_iter=2,000. StandardScaler is fitted on discovery only. Signature classification uses the four signed expression-profile scores; pathway/null models use the same logistic recipe.

| model | accuracy | balanced_accuracy | macro_f1 | CMS1_recall | CMS2_recall | CMS3_recall | CMS4_recall |
| --- | --- | --- | --- | --- | --- | --- | --- |
| majority | 0.4273 | 0.2500 | 0.1497 | 0.0000 | 1.0000 | 0.0000 | 0.0000 |
| batch_only | 0.3636 | 0.2746 | 0.2697 | 0.1333 | 0.6596 | 0.1765 | 0.1290 |
| expression_logistic | 0.9182 | 0.9169 | 0.9314 | 1.0000 | 0.9787 | 0.8824 | 0.8065 |
| expression_nearest_centroid | 0.8818 | 0.8936 | 0.8898 | 0.9333 | 0.8936 | 0.9412 | 0.8065 |
| signed_signature_scores | 0.8909 | 0.9252 | 0.9021 | 1.0000 | 0.8298 | 1.0000 | 0.8710 |
| pathway_logistic | 0.9455 | 0.9516 | 0.9582 | 1.0000 | 1.0000 | 1.0000 | 0.8065 |
| shuffled_pathway_17 | 0.9091 | 0.9049 | 0.9069 | 1.0000 | 0.9574 | 0.8235 | 0.8387 |
| shuffled_pathway_29 | 0.9364 | 0.9397 | 0.9400 | 1.0000 | 0.9787 | 0.9412 | 0.8387 |
| shuffled_pathway_43 | 0.9000 | 0.9035 | 0.9143 | 1.0000 | 0.9574 | 0.8824 | 0.7742 |

The three shuffled-membership controls test whether curated pathway assignments improve reproduction beyond arbitrary expression pooling at the same feature count. Nulls may retain gene signal, and overlapping pathways make this a restricted control. A high subtype score is **not novel subtype discovery**: CMS labels were themselves derived from expression, and this cohort contributed to the published CMS work. The original discovery/validation grouping is useful, but does not make the CMS labels independently assigned from non-expression evidence. A new external cohort is required for stronger claims. Accuracy deltas against expression logistic are described in `paired_bootstrap_accuracy_delta.csv`; within-class patient resampling on one panel does not establish external generalization or undo label circularity.

## Direct drug-target connectivity of signatures

| subtype | selected_genes | selected_genes_with_targets | unique_connected_drugs | interpretation |
| --- | --- | --- | --- | --- |
| CMS1 | 100 | 21 | 74 | connectivity only; no action direction, vulnerability or efficacy inferred |
| CMS2 | 100 | 18 | 117 | connectivity only; no action direction, vulnerability or efficacy inferred |
| CMS3 | 100 | 21 | 65 | connectivity only; no action direction, vulnerability or efficacy inferred |
| CMS4 | 100 | 10 | 81 | connectivity only; no action direction, vulnerability or efficacy inferred |

Only PrimeKG `drug_protein` display role `target` is counted; reciprocal rows and aliases are deduplicated by stable IDs. These are connectivity counts, **not candidate-drug rankings**. Increased/decreased expression does not establish tumor dependency, required pharmacological direction, or therapeutic benefit. A final drug–subtype predictor still needs suitable supervision and independent response evaluation.

## What this supports and does not support

The appropriate conclusion is whether these CMS expression profiles are a reproducible, identifier-audited molecular input. Compare the measured profile consistency and classifier/null results before deciding whether pathway structure adds value. Better CMS reproduction alone cannot establish drug prediction, GNN superiority, genotype causality, patient survival benefit, or clinical efficacy. It also cannot establish a new contribution beyond existing CMS classification.

This cohort has selected mutation-status fields, not the full mutation matrix used in ECoGraph. CMS1–4 must remain separate from GI.CIN/GI.HM-indel/GI.HM-SNV/GI.GS; there is no assumed one-to-one crosswalk. Bulk expression contains immune/stromal composition, particularly relevant to mesenchymal/immune signatures. Batch, stage and location can confound effects despite adjustment; OLS is a descriptive model with limited validation group sizes. FDR and bootstrap intervals do not correct provenance or selection bias. The filename suggests fRMA preprocessing, but raw CEL-level preprocessing was not independently reproduced. Four clinical records have stage0; their interpretation warrants upstream review and no assumption that these are matched normal controls is made.

## Next decision

If profiles replicate, retain them as a **new CMS representation benchmark**, not a drop-in proof of better drug inference. A future GNN can integrate these profiles through mapped genes, but should first demonstrate a benefit over expression/pathway baselines on a clearly defined held-out outcome. For subtype–drug links, obtain independent cell-model/experimental response evidence with defensible CMS assignments. Graph indication learning and literature review alone do not verify therapeutic benefit; literature retrieval should audit support, contradiction, and prior knowledge with citations.

## Files and reproduction

Run `python3 analysis/crcsc_experiment.py`; `--output` must not already exist. Dependencies: pandas, NumPy, SciPy, scikit-learn, threadpoolctl. Official annotation is a local prerequisite; no runtime download occurs. Source files are unchanged (SHA-256 audit: **True**). See `manifest.json` for versions and hashes.

Outputs: `sample_audit.csv`, `cohort_counts.csv`, `probe_mapping_audit.csv.gz`, `gene_mapping_audit.csv`, `gene_profiles.csv.gz`, `selected_gene_signatures.csv`, `pathway_profiles.csv.gz`, `pathway_inventory.csv`, `profile_replication.csv`, `signature_overlap.csv`, `signature_drug_connectivity.csv`, `classification_metrics.csv`, `validation_predictions.csv`, `confusion_*.csv`, `batch_by_split.csv`, `paired_bootstrap_accuracy_delta.csv`, `manifest.json`, `findings.md`.

Primary context: [CMS study](https://www.nature.com/articles/nm.3967), [GSE39582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582). Local file integrity is checked; exact download provenance of the three supplied files and upstream preprocessing remains to be documented by the team.
