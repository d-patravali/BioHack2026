# CRC GNN proof-of-concept: findings

**Verdict:** the pipeline generates traceable association hypotheses and supports quantitative sensitivity tests. Drug popularity exceeds every fitted GNN in this run. The tested profiles do not fix the shared top candidate group. No subtype response accuracy or clinical efficacy was evaluated.

## What was run

The existing notebook architecture, width 48, optimizer and 40-epoch budget were retained. 8 fits cover the legacy model (seed42), direct-target-only pooling and a retrained drug-only molecular-context ablation (initialization seeds [42, 43, 44]), and no-PPI (seed42). All use exactly the same disease split, positive examples, epoch-specific sampled unknown training pairs and validation/test pairs. Legacy AP matches the historical notebook to displayed precision: **True**. No original notebooks or outputs changed; all preservation checks passed.

There are 6,233 eligible indications: 4,229 train, 1,151 validation and 853 test; 633/136/136 diseases respectively. Test labels comprise 853 positives and 853 unknowns. Drug popularity uses training counts only. The fixed evaluation retains all 6,282 original drugs, including 592 with no direct target message in corrected models; a common target-supported pair-subset metric is also saved. Hypothesis rankings are restricted to 5,690 direct-target-supported drugs. No therapeutic edges enter molecular adjacency. Source/family overlap and false negatives remain possible.

## 1. Demonstrated predictive performance

| condition | runs | mean_test_AP | min_test_AP | max_test_AP | mean_disease_macro_AP |
| --- | --- | --- | --- | --- | --- |
| drug_popularity | 1 | 0.8941 | 0.8941 | 0.8941 | 0.9218 |
| legacy_all_roles | 1 | 0.8190 | 0.8190 | 0.8190 | 0.8580 |
| target_only | 3 | 0.8254 | 0.8158 | 0.8337 | 0.8490 |
| drug_only | 3 | 0.8257 | 0.8202 | 0.8314 | 0.8442 |
| no_ppi | 1 | 0.8125 | 0.8125 | 0.8125 | 0.8312 |

Drug popularity exceeds every fitted GNN in this run. Across initialization seeds, target-only minus drug-only test AP is **-0.0003 on average**, ranging **-0.0043 to +0.0023**. This comparison measures whether removing molecular disease context degrades prediction. The single-seed PPI comparison is exploratory, not evidence of reliable graph benefit. These ranges vary initialization on one split, not independent datasets or confidence intervals. Every selected checkpoint reached the epoch-budget boundary: **True**. A fixed budget does not rule out improvement with other training.

## 2. Why subtype rankings are nearly identical

Two effects are directly observed:

- **Drug-level variation dominates.** Across-drug main effects account for the following percentage of total score variance. The rest includes subtype-wide offsets and drug-by-subtype interactions; this decomposition is descriptive, not a causal attribution.

| profile | mean | min | max |
| --- | --- | --- | --- |
| frequency | 99.9995 | 99.9989 | 99.9999 |
| selective_frequency | 99.9679 | 99.9431 | 99.9826 |
| uniform_union | 99.9998 | 99.9995 | 100.0000 |

- **Identical target neighborhoods are structurally indistinguishable to this model.** Observed top-score tie-group sizes across target-only seeds and actual profiles are **[29]**; inclusive top20 counts are **[29]**. All tie-inclusive selective-profile candidates have NR3C1 (NCBI2908) as their only recorded direct target: **True**. That gene has one PrimeKG node. Identical target neighborhoods necessarily produce identical drug embeddings and scores in this architecture, which has no compound-specific features. These are recorded connections, not evidence that modulating NR3C1 benefits CRC. `tie_aware_top_candidates.csv` reports minimum tied rank and full group size rather than inventing an ordering among those compounds.

The broad union profiles and mean pooling also discard mutation frequency and average many overlapping genes. This is consistent with the collapse, while the controls quantify the limited dependence rather than proving a single cause.

## 3. Did weighted or selective profiles help?

Profiles are all-observed uniform weights; normalized observation frequencies; and a prespecified selective rule (observed in >=2 patients and >=10% of the subtype, capped at top200 by frequency). They use no indication test labels or favored drugs. Selection is descriptive, not burden-adjusted enrichment.

| subtype | uniform_union | frequency | selective_frequency |
| --- | --- | --- | --- |
| subtype:GI.CIN | 9763 | 9763 | 22 |
| subtype:GI.GS | 3741 | 3741 | 57 |
| subtype:GI.HM-SNV | 11421 | 11421 | 200 |
| subtype:GI.HM-indel | 15580 | 15580 | 200 |

Across six subtype pairs and three target-only seeds:

| profile | mean_top20_overlap | mean_score_Spearman | mean_pairwise_score_RMS |
| --- | --- | --- | --- |
| uniform_union | 20.0000 | 1.0000 | 0.0062 |
| frequency | 20.0000 | 1.0000 | 0.0109 |
| selective_frequency | 20.0000 | 0.9999 | 0.0955 |
| permuted_gene_control | 20.0000 | 0.9999 | 0.1359 |

All raw top20 subtype-pair overlaps equal20: **True**. The tested profiles do not fix the shared top candidate group. Numerical score changes are a sensitivity result, not validated improvement. `profile_summary.csv` audits cutoff ties: the selective HMSNV cap in this run cuts through a 412-gene frequency tie. Its deterministic tie choice is especially unstable with only six patients. Frequencies are not driver probabilities and the procedure does not resolve mutation burden, gene length, passenger flags or missing variant consequences.

## 4. Ablations and interpretation

- Identical profiles and zero molecular profiles give identical subtype predictions, as expected.
- Swapping profiles swaps their prediction vectors exactly, confirming computational dependence on supplied context rather than the subtype name.
- Permuting gene identities changes scores and can increase numerical separation while preserving the same top drug group. Arbitrary inputs can differentiate outputs; differentiation alone is not biological validation.
- The retrained drug-only ablation matches the context model closely on held-out indication AP. This undermines a claim that the current supervised task learns useful disease-specific molecular context.

Subtype-only weighting is applied after training: subtypes have no indication labels, and gene->disease pooling has no feedback into gene/drug embeddings. Changing their aggregation rows cannot change the disease evaluation AP. We therefore do not claim weighting improved supervised performance. Profile effects are a separate inference-sensitivity experiment.

## 5. Responsible proof-of-concept claim and next step

**Demonstrated:** reproducible compact graph inputs; original-model reproduction; corrected target semantics; quantitative profile perturbations; exact-tie auditing; fair baseline/ablation comparisons. **Not demonstrated:** superior indication prediction, useful subtype-specific top candidates, subtype response accuracy, clinical efficacy, or literature novelty.

The defensible demo is an investigative agent exposing hypotheses **alongside** their tied target groups, profile evidence, and failed baseline comparison. Do not market these outputs as validated or novel treatments. Prioritize an evaluation that rewards disease/subtype discrimination over shared drug popularity, and audit the contrast diagnostic with independent evidence before further tuning. A more distinct list is not itself an improvement. Existing disease-level validation cannot establish which of these four unlabeled subtype hypotheses is correct.

Raw logits and within-drug subtype contrasts in `hypothesis_rankings.csv` are explicitly unvalidated. The contrasts remove a shared drug preference for inspection; they are not a second validated predictor and should not be selected because they look more interesting. Control-profile rankings are diagnostics, not hypotheses to recommend. Graph nonmembership does not establish novelty in the literature. `recorded_indication_for_crc_lexical_candidate` only audits this graph's broad name-search candidates, not curated subtype efficacy.

## Reproduction

`python3 analysis/crc_experiments.py` creates a new timestamped run. Requires the local HGNC TSV supplied with `--hgnc` (default `data/external/hgnc/hgnc_complete_set.txt`); the experiment runner does not download it. `python3 analysis/crc_report.py reports/crc_experiments/20261108_034016` regenerates this summary from completed results without retraining. See `crc_experiment_report.md` for extended methods and `provenance.json` for hashes, versions and preservation evidence.
