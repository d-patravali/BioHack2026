# CRC hypothesis-generation experiments

## Scope and minimal changes

Separate scripted experiments reuse the notebook's architecture (width 48, one PPI propagation step, drug/disease gene pooling and the same decoder), 40-epoch default, AdamW settings and disease-held-out protocol. No original notebooks or saved outputs were modified. No full merged graph artifact was available locally, so compact model inputs were reconstructed from the unchanged raw PrimeKG and ECoGraph matrix using the notebook's HGNC mapping rule. The HGNC reference is a new downloaded snapshot, fingerprinted in provenance; it need not equal the teammate's historical snapshot.

The changes are weighted/selective **inference-time** subtype profiles, direct-target filtering as a separately identified model condition, and context/PPI ablations. Subtypes have no supervised indication examples. In this feed-forward model, changing only their outgoing aggregation rows cannot alter training of labeled diseases, so retraining each profile would add no experimental information. The existing model is retained rather than replaced.

## Inputs and fixed evaluation

Input counts: `{"raw_rows": 8100498, "genes": 27671, "drugs_all_roles": 6282, "drugs_direct_target": 5690, "zero_target_drugs_retained_for_matched_evaluation": 592, "mapped_matrix_genes": 18575, "known_indication_pairs": 9388, "eligible_positive_pairs": 6233, "direct_target_pairs": 16380, "all_role_drug_gene_pairs": 25468, "typed_patients": 339}`.

Disease split (seed 42): `{"train": {"diseases": 633, "positive_pairs": 4229}, "validation": {"diseases": 136, "positive_pairs": 1151}, "test": {"diseases": 136, "positive_pairs": 853}}`. Every condition uses the exact same train positives, per-epoch sampled training pairs, validation pairs and test pairs. Model seeds [42, 43, 44] vary initialization, not split or patient sampling. Validation AP alone chooses an epoch. Legacy and no-PPI are single-seed diagnostics. Original target/reverse rows are canonicalized; indication/off-label/contraindication labels are excluded from molecular adjacency and from sampled unknown pools. Unknown pairs are not biological negatives. Only exact node identities are split; source overlap and disease-family equivalence remain unresolved.

The candidate universe for disease evaluation deliberately retains **all legacy drug-protein-connected drugs** so correcting target semantics does not silently change the task. Target-only models therefore have 592 drugs with zero target messages; `test_target_supported_ap` separately reports the same pair subset supported by actual target edges. CRC hypothesis rankings use only 5690 direct-target-supported drugs. AP is computed on a 1:1 positive/unknown sample, not the full-catalog retrieval task or therapeutic efficacy. Drug popularity is computed on training indications only.

## Demonstrated indication-recovery performance

| condition | runs | mean_test_ap | min_test_ap | max_test_ap | disease_macro_ap |
| --- | --- | --- | --- | --- | --- |
| drug_popularity | 1 | 0.8941 | 0.8941 | 0.8941 | 0.9218 |
| legacy_all_roles | 1 | 0.8190 | 0.8190 | 0.8190 | 0.8580 |
| target_only | 3 | 0.8254 | 0.8158 | 0.8337 | 0.8490 |
| drug_only | 3 | 0.8257 | 0.8202 | 0.8314 | 0.8442 |
| no_ppi | 1 | 0.8125 | 0.8125 | 0.8125 | 0.8312 |

These numbers are the measured results; a changed subtype ranking cannot compensate for poor held-out performance. Initialization ranges on one fixed split are not confidence intervals or external validation. Disease-macro AP reduces domination by heavily annotated diseases but is still a sampled annotation-recovery metric. Historical notebook values (GNN AP 0.819, popularity AP 0.894) are context only; reconstructed inputs/software can change exact values.

## Why the original subtype rankings can collapse

The original code ignores `mutation_fraction` and averages broad, overlapping observed-gene unions. It also collapses disease-gene associations and subtype mutation observations into the same aggregation rule. Gene averaging, PPI smoothing, a drug-driven decoder and annotation popularity can all reduce subtype dependence; the experiments below measure their consequences rather than assuming a cause. Both the original and target-only drug pooling still omit action direction and variant function.

Profile construction is prespecified, with no disease-test or desired-drug tuning: (1) every observed mapped gene equally weighted; (2) all positive genes weighted by subtype observation frequency; (3) retain genes observed in at least two patients and at least 10% of the subtype, cap at 200 highest frequencies, and normalize frequency weights. This is descriptive selection, **not burden-adjusted enrichment**. A fixed node-index tie rule is used; cutoff ties are audited. Six-patient GI.HM-SNV remains uncertain. `profile_summary.csv` reports actual gene counts and effective weighted sizes.

## Subtype sensitivity: differentiation is not validation

The table averages the six subtype pairs and available initialization seeds for the target-only model. Spearman uses tied scores; deterministic node-index tie breaking is used for top20 membership, so rank-list differences near ties can be arbitrary.

| profile | mean_top20_overlap | mean_score_spearman | mean_rms_difference |
| --- | --- | --- | --- |
| uniform_union | 20.0000 | 1.0000 | 0.0062 |
| frequency | 20.0000 | 1.0000 | 0.0109 |
| selective_frequency | 20.0000 | 0.9999 | 0.0955 |
| identical_pooled_control | 20.0000 | 1.0000 | 0.0000 |
| zero_context_control | 20.0000 | 1.0000 | 0.0000 |
| permuted_gene_control | 20.0000 | 0.9999 | 0.1359 |
| swapped_subtype_control | 20.0000 | 0.9999 | 0.0955 |

For the following table, each condition is compared with that same subtype's selective-frequency profile using identical model weights:

| comparison_profile | mean_top20_overlap | mean_score_spearman | mean_rms_change |
| --- | --- | --- | --- |
| uniform_union | 20.0000 | 1.0000 | 0.0648 |
| frequency | 20.0000 | 1.0000 | 0.0606 |
| identical_pooled_control | 20.0000 | 1.0000 | 0.0567 |
| zero_context_control | 20.0000 | 0.9983 | 1.4203 |
| permuted_gene_control | 20.0000 | 0.9998 | 0.1558 |

Controls: identical pooled profiles must produce identical scores; zero profiles remove subtype molecular input; a fixed gene permutation preserves each profile's weights/size but destroys gene identity correspondence; swapping the four profiles must swap predictions exactly. Assertions verify identical-profile and swap behavior. The retrained drug-only condition removes all disease molecular context, testing whether labeled disease gene neighborhoods improve the supervised task. The no-PPI condition is a single-seed graph-structure diagnostic. These are complementary tests: perturbation sensitivity shows computational dependence, whereas predictive improvement over controls would support task-relevant information. Random profiles can also generate different rankings; differentiation alone is not success.

`embedding_diagnostics.csv` compares subtype embedding cosines and between-subtype score RMS relative to across-drug score variation. `scores_*.npz` retain full candidate scores for auditing. `hypothesis_rankings.csv` provides raw logits and a separately labeled within-drug subtype contrast (score minus that drug's mean over four subtypes). Contrast is a diagnostic to remove a shared drug preference, **not a second validated predictor**. Hypotheses are not calibrated probabilities, clinical recommendations, or evidence of differential response. No independent subtype response labels or prospective novelty validation are available.

## Interpretation and next steps

Use the measured performance and controls together. If weighting/selectivity changes rankings but not biological predictive evidence, report improved controllability/sensitivity only. If drug-only or popularity matches/exceeds the GNN, acknowledge that the present experiment has not demonstrated a GNN advantage. If permuted profiles produce comparable differentiation, do not infer that differentiated rankings are biologically meaningful. A no-PPI comparison at one seed is exploratory and cannot establish a reliable graph benefit.

The defensible proof of concept is traceable molecular input -> graph model -> explicitly unvalidated association hypotheses, with honest baselines and sensitivity audits. Before selecting hypotheses for further investigation, inspect variant consequences, driver/passenger context, drug action, known CRC indications, assay evidence and source independence. Disease-level evaluation does not validate the four unlabeled subtype predictions. More distinct lists are not a scientific objective on their own.

## Reproduction and preservation

`python3 analysis/crc_experiments.py` creates a new timestamped directory; `--output` must name a nonexistent directory. Defaults are 40 epochs and model seeds 42,43,44 on CPU with two threads and deterministic algorithms. The script never executes notebook cells, downloads data itself, changes original outputs, or writes model results into the original prediction directories. All protected inputs, notebooks and preexisting report/artifact files were fingerprinted before and after and remained unchanged. Versions, hashes, split IDs, sampled pairs, training curves and diagnostic scores are saved. Dependencies are existing pandas, NumPy, SciPy, scikit-learn and PyTorch.
