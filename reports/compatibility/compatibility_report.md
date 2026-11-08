# ECoGraph–PrimeKG compatibility assessment

## Verdict: VIABLE WITH MODIFICATIONS

**Technical gene-level integration is plausible, but the supplied Neptune export is not ready to merge.** Its sample–gene file is a two-byte blank line with no header or observations. Subtype molecular networks cannot be computed from those eight exports alone. All subtype gene-set numbers below come from a separately verified, local binary mutation matrix, not from the empty edge file. No source files were repaired, no graph was merged, and no GNN, agent, drug ranking, or efficacy prediction was implemented.

The conditional verdict means: repair/validate the sample–gene export, remove a metadata-as-gene artifact, preserve unharmonized pathways separately, and define a defensible prediction/evaluation target. **Predictive feasibility and therapeutic efficacy remain unestablished.**

## Research question and inputs

Can GI.CIN, GI.HM-indel, GI.HM-SNV, and GI.GS molecular observations be linked to PrimeKG for downstream drug–subtype hypothesis generation? This assessment separates technical compatibility, molecular associations, target connectivity, therapeutic supervision, and efficacy. Inputs are `data/raw/primekg/kg.csv` and `ECoGraph/data/neptune_csv/` relative to the project root. `source_inspection.txt` records headers and example records; `source_manifest.json` records sizes and SHA-256 digests.

### ECoGraph export structure and integrity

| file | rows | status | duplicate_ids | duplicate_rows |
| --- | --- | --- | --- | --- |
| vertices_gene.csv | 21333 | parsed | 0 | 0 |
| vertices_subtype.csv | 4 | parsed | 0 | 0 |
| vertices_sample.csv | 404 | parsed | 0 | 0 |
| vertices_age.csv | 2 | parsed | 0 | 0 |
| edges_has_gene_part001.csv | 0 | empty_no_schema | 0 | 0 |
| edges_has_subtype.csv | 339 | parsed | 0 | 0 |
| edges_in_pathway.csv | 35818 | parsed | 0 | 0 |
| edges_is_of_age.csv | 404 | parsed | 0 | 0 |

Declared node labels are Gene, Sample, Subtype, and Age. Gene vertices contain `~id`, `~label`, `symbol:String`, `exp_hr_coeff:String`, and `hr_pval:String`. Symbols, not Entrez/Ensembl IDs, define `gene:SYMBOL` identifiers. Sample vertices use `sample:TCGA-...` patient-level IDs plus a patient barcode and sex. Subtypes are `subtype:GI.*`; age IDs are `age:<=50` and `age:>50`. The two survival fields contain literal missing markers, counted separately from blanks in `ecograph_integrity.csv`; they are not subtype evidence or drug-response labels.

Edges use `~id`, `~from`, `~to`, and `~label`. Verified directions are Sample → Subtype (`HAS_SUBTYPE`), Gene → Pathway (`IN_PATHWAY`), and Sample → Age (`IS_OF_AGE`). The sample–gene file has no schema, so its direction and semantics cannot be verified from the export. The accompanying notebook intends Sample → Gene (`HAS_GENE`) with nonzero mutation observations, but this is code intent, not populated edge evidence.

The 339 subtype edges give one assignment each to 224 GI.CIN, 60 GI.HM-indel, 6 GI.HM-SNV, 49 GI.GS samples. There are 404 sample vertices; 65 have no exported subtype. All existing subtype and age edge endpoints resolve to supplied vertices. All pathway source genes resolve; **every pathway destination is undeclared** because no `vertices_pathway.csv` is supplied: 35,818 edges reference 2,036 distinct pathway endpoints. None of the nonempty ECoGraph files has duplicated IDs or rows or blank cells. Missing marker tokens still occur in gene properties.

**Metadata artifact:** gene vertices contain `gene:patient_id`. The mutation matrix has a metadata column `patient_id`, and the notebook's `gene_columns` exclusion list omits it. This is not a biological gene and is excluded from all biological denominators and subtype sets, while retained in the mapping audit. There are 21,333 exported gene vertices and **21,332 eligible gene symbols** (including symbols unobserved in a particular subtype).

### Supplementary recovery, not silent export repair

`ECoGraph/data/coad_mut_data.csv` contains 404 distinct patient records and 21332 binary gene columns; patient IDs equal the first 12 characters of tumor barcodes. Its patient IDs match all 404 sample vertices. After excluding `Tumor_Sample_Barcode` and `patient_id`, gene columns match eligible exported symbols exactly. Values were verified to be only 0/1. Every exported subtype pair agrees with `coad_subtype_data.csv`; two additional GI.CIN assignments in that source lack corresponding exported samples and are not introduced.

The supplementary analysis joins the binary matrix to the **exported** 339 subtype assignments, counts each patient–gene pair once, and interprets 1 as an observed binary mutation flag according to the local processing documentation. It does not count variants, infer activating/inactivating status, or regard a 0 as proof that the gene was comprehensively tested. 151776 observed pairs occur across all 404 matrix patients. No reconstruction is written into the source export. Export-only subtype molecular metrics remain **NA**, not zero, in `neptune_only_subtype_status.csv`.

## PrimeKG schema, identifiers, and relevant relationships

The raw CSV scan finds 8,100,498 ordered rows, 129,375 node indices, 10 node types and 30 relation labels. Node types are anatomy, biological_process, cellular_component, disease, drug, effect/phenotype, exposure, gene/protein, molecular_function, pathway. Each edge row has `relation`, `display_relation` and endpoint `index`, `id`, `type`, `name`, `source` fields. All identifier fields are read as strings. No blank cells or conflicting metadata per node index were found. There are 126 node rows sharing `(type, source, id)` (122 gene rows, 4 anatomy rows); exact-ID aliases are counted once for biological matching, without modifying original nodes.

Gene/protein nodes use NCBI IDs (27,671 original nodes, 27,610 distinct IDs). Pathways use Reactome stable `R-HSA-...` identifiers (2,516). Drugs use DrugBank `DB...` IDs (7,957). Diseases use MONDO and MONDO_grouped IDs (17,080 nodes). `primekg_node_types.csv` and `primekg_relations.csv` contain the complete observed namespace and typed relation inventory.

Crucially, `drug_protein` has four display labels: **target, enzyme, carrier, transporter**. Only `display_relation == target` is counted as direct drug–target connectivity. Enzymes, carriers, and transporters are excluded from target counts. Drug→Gene and Gene→Drug rows are reciprocal representations, verified and deduplicated. Likewise, `indication` rows are deduplicated by DrugBank drug and disease ID; `contraindication` and `off-label use` are kept distinct and are not therapeutic positives in these counts.

| display_relation | unique_drug_gene_pairs |
| --- | --- |
| carrier | 864 |
| enzyme | 5317 |
| target | 16380 |
| transporter | 3092 |

There are **16,380 direct target pairs**, involving **2,836 genes and 5,690 drugs**. Context relations also include `protein_protein`, `disease_protein`, `pathway_protein`, GO functional annotations, anatomy expression, and positive/negative disease–phenotype links. These offer possible message-passing context, not proof of drug response.

## Gene identifier compatibility

ECoGraph lacks stable gene identifiers. Mapping uses exact, case-sensitive gene symbols against PrimeKG gene/protein `name`, accepting a match only if it resolves to one NCBI ID. Only leading/trailing whitespace may be normalized, with a separate status. We do not uppercase everything, remove punctuation, use fuzzy names, or resolve aliases automatically. Multiple indices with the same NCBI ID are one biological candidate; a symbol mapping to multiple NCBI IDs is ambiguous and is not accepted. Accepted symbol matches are provisional identifier bridges, not an authoritative HGNC crosswalk. A curated stable-ID crosswalk is recommended before integration.

| raw_export_gene_vertices | eligible_ecograph_genes | primekg_gene_nodes | unique_primekg_genes | matching_ecograph_genes | matching_unique_ncbi_ids | ecograph_mapping_pct | raw_export_mapping_pct | primekg_in_ecograph_pct | ambiguous_ecograph_genes | unmatched_eligible_genes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 21333 | 21332 | 27671 | 27610 | 18909 | 18909 | 88.64 | 88.64 | 68.49 | 0 | 2423 |

The main percentage excludes the known `patient_id` artifact; the CSV also reports raw-export coverage so the denominator is explicit. `identifier_mapping_audit.csv` and `gene_mapping_audit.csv` retain exact, normalized, ambiguous, unmatched, and excluded rows, accepted IDs, and all candidate PrimeKG indices. A high overlap here establishes useful bridges, not subtype-specific therapeutic relevance.

## Subtype molecular associations and mapping

**All tables in this section are supplementary matrix-derived.** “Associated” means observed in at least one assigned sample, not enriched, causal, a driver, or exclusive biology. Sample–gene observations are unique binary patient–gene pairs.

| subtype | samples | associated_genes | sample_gene_observations | mapped_genes | mapping_pct | unmatched_genes | ambiguous_genes | exclusive_observed_genes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GI.CIN | 224 | 10590 | 24376 | 9909 | 93.57 | 681 | 0 | 360 |
| GI.HM-indel | 60 | 17003 | 76068 | 15829 | 93.10 | 1174 | 0 | 2966 |
| GI.HM-SNV | 6 | 12333 | 25059 | 11603 | 94.08 | 730 | 0 | 603 |
| GI.GS | 49 | 4022 | 5565 | 3779 | 93.96 | 243 | 0 | 75 |

Per-gene counts and frequencies are provided in the four `subtype_gene_sets/*.csv` files. Frequency is the number of subtype patients with an observed gene flag divided by subtype sample count. Set size and exclusivity are highly sensitive to both sample count and mutation burden.

Selected concrete gene examples (descriptive frequencies, not drug-response biomarkers):

| symbol | NCBI ID | GI.CIN observed % | GI.HM-indel observed % | GI.HM-SNV observed % | GI.GS observed % | PrimeKG direct target |
| --- | --- | --- | --- | --- | --- | --- |
| APC | 324 | 82.59 | 38.33 | 100.00 | 81.63 | False |
| TP53 | 7157 | 75.89 | 28.33 | 33.33 | 14.29 | True |
| KRAS | 3845 | 43.30 | 28.33 | 33.33 | 71.43 | True |
| BRAF | 673 | 4.91 | 60.00 | 33.33 | 10.20 | True |
| PIK3CA | 5290 | 25.00 | 40.00 | 100.00 | 44.90 | True |

For example, APC can map into PrimeKG even though no direct drug-target edge is recorded for it; a recorded target connection for TP53 or KRAS does not establish that any drug is appropriate for the observed variant. Missing target edges are also not proof that a gene is untargetable.

### Mutation burden and enrichment screens

| subtype | samples | minimum | q25 | median | q75 | maximum | mean |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GI.CIN | 224 | 45 | 82.75 | 101.50 | 130.25 | 267 | 108.82 |
| GI.HM-indel | 60 | 357 | 928.50 | 1112.50 | 1605.50 | 3326 | 1267.80 |
| GI.HM-SNV | 6 | 2572 | 2933.00 | 4285.00 | 5116.50 | 6041 | 4176.50 |
| GI.GS | 49 | 53 | 89.00 | 106.00 | 138.00 | 225 | 113.57 |

These burdens count distinct flagged genes per sample, **not mutations per megabase**. Hypermutated groups have more opportunities to accumulate passenger flags; GI.HM-SNV has only six assigned samples. We therefore do not treat a union of observed genes as a subtype-specific molecular signature.

We executed two explicitly exploratory one-sided tests for increased binary prevalence, comparing each subtype with all other typed patients: exact Fisher/hypergeometric tests, and a burden-stratified Cochran–Mantel–Haenszel normal-approximation screen over pooled burden quintiles. The latter sums observed-minus-expected counts and conditional variances in strata containing both groups; strata without overlap contribute no information. We suppress the burden screen unless each group has at least 20 samples and at least half of each group's samples fall in mixed strata. These are minimal operational guardrails, not proof of test validity. Disabled tests receive p=1 before BH; their summary results are NA. Unfiltered approximate p-values are retained only for diagnostics. BH correction is applied across all 4 × 21,332 hypotheses separately for each test family. The joint screen requires both q-values <0.05 and greater observed subtype prevalence. All tests and support counts are exported; this is not a validated enrichment procedure.

| subtype | raw_fisher_bh_hits | burden_screen_eligible | subtype_samples_in_mixed_strata | other_samples_in_mixed_strata | burden_cmh_bh_hits | joint_exploratory_screen_candidates | screen_candidates_mapped | screen_candidates_mapped_pct | screen_candidates_with_direct_targets | validated_enriched_genes | interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GI.CIN | 2 | True | 224 | 115 | 8.00 | 1.00 | 1.00 | 100.00 | 1.00 | NA | No validated enriched set: asymptotic burden screen, sparse groups, no external replication |
| GI.HM-indel | 6353 | False | 60 | 8 | NA | NA | NA | NA | NA | NA | No validated enriched set: asymptotic burden screen, sparse groups, no external replication |
| GI.HM-SNV | 5247 | False | 6 | 62 | NA | NA | NA | NA | NA | NA | No validated enriched set: asymptotic burden screen, sparse groups, no external replication |
| GI.GS | 3 | True | 49 | 222 | 0.00 | 0.00 | 0.00 | NA | 0.00 | NA | No validated enriched set: asymptotic burden screen, sparse groups, no external replication |

The HM-indel burden strata contain only eight comparison patients, so that burden screen is disabled; GI.HM-SNV also fails the sample-size and comparison-support guardrails. Naive Fisher testing produces thousands of hypermutated-group hits, illustrating why unadjusted prevalence must not be presented as subtype-specific biology. The CMH approximation can still be poor for sparse genes; burden quintiles provide incomplete adjustment and include each tested gene in the total burden. No callable-region, gene-length, variant-consequence, sequencing-depth, treatment, or external replication data are used. Thus any eligible joint candidates are **exploratory associations, not validated subtype-enriched genes**. Required validated-enrichment counts remain NA. Their mapped/target counts are separately labeled exploratory, rather than promoted to reliable enriched-gene coverage.

## Pathway compatibility

ECoGraph's 2,036 endpoint IDs are collection-prefixed names: 1,567 REACTOME, 284 BIOCARTA, and 185 KEGG. Example: `pathway:BIOCARTA_41BB_PATHWAY`. PrimeKG has 2,516 Reactome stable IDs such as `R-HSA-109581` (Apoptosis). ECoGraph does not supply corresponding stable database IDs or a name-to-ID crosswalk.

| ecograph_pathways | primekg_pathways | confident_matches | confident_mapping_pct | unmatched_without_confident_mapping | ambiguous_label_candidates | unique_unverified_label_candidates | biological_overlap_pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2036 | 2516 | 0 | 0.00 | 2036 | 4 | 1542 | NA |

No pathway is accepted solely on a similar or normalized name. The audit suggests candidate Reactome label correspondences by stripping `REACTOME_` and normalizing punctuation/case; these are explicitly unverified and never used to merge entities or claim overlap. **Zero confidently established mappings is not evidence of zero biological overlap. Actual pathway overlap percentage is not reliably computable from supplied identifiers.** BIOCARTA and KEGG pathways must remain separate from Reactome entities unless a defensible equivalence relation is supplied. A future crosswalk should preserve source collection, version, stable ID, and mapping provenance.

## Direct drug-target connectivity by subtype

| subtype | mapped_genes | mapped_genes_with_direct_targets | direct_target_pct_of_mapped_genes | unique_connected_drugs | unique_direct_drug_gene_pairs |
| --- | --- | --- | --- | --- | --- |
| GI.CIN | 9909 | 1555 | 15.69 | 4231 | 10657 |
| GI.HM-indel | 15829 | 2461 | 15.55 | 5373 | 15163 |
| GI.HM-SNV | 11603 | 1819 | 15.68 | 4667 | 12111 |
| GI.GS | 3779 | 601 | 15.90 | 2081 | 4395 |

These are direct `drug_protein`/`target` edges, not indirect graph paths. Drug counts are distinct DrugBank IDs; gene coverage counts distinct ECoGraph symbols resolving to unambiguous PrimeKG IDs. The CSV separately records distinct mapped NCBI IDs. A connected drug may target an observed gene without modulating the relevant variant, acting in the needed direction, reaching the tissue, or improving disease. Connectivity, indication, prediction, and actual efficacy are separate quantities. There are no observed drug–GI subtype therapeutic labels in these inputs.

## Cross-subtype comparison

| subtype_a | subtype_b | genes_a | genes_b | overlap | jaccard |
| --- | --- | --- | --- | --- | --- |
| GI.CIN | GI.HM-indel | 10590 | 17003 | 9829 | 0.55 |
| GI.CIN | GI.HM-SNV | 10590 | 12333 | 7791 | 0.51 |
| GI.CIN | GI.GS | 10590 | 4022 | 2894 | 0.25 |
| GI.HM-indel | GI.HM-SNV | 17003 | 12333 | 11306 | 0.63 |
| GI.HM-indel | GI.GS | 17003 | 4022 | 3831 | 0.22 |
| GI.HM-SNV | GI.GS | 12333 | 4022 | 3190 | 0.24 |

**2,414 observed genes are shared across all four subtypes.** Exclusives are listed by count above; gene-set CSVs allow identification. `subtype_exploratory_screen_overlap.csv` repeats comparisons for joint screen candidates, without treating them as validated enrichment. Different observed unions and frequencies can give different input representations, but may primarily encode mutation burden and sample size. They do not establish different therapeutic vulnerabilities. Weighted prevalence representations and burden-controlled validation should precede any claim of meaningful subtype drug-ranking differences. Mapping percentages and direct-target percentages also differ, as quantified above.

## Therapeutic supervision and colorectal disease context

PrimeKG contains **9,388 unique known indication pairs**, spanning **1,801 drugs and 1,363 diseases**. These provide a technically usable supervised drug–disease relation, although abundance alone does not establish adequacy for the specific CRC subtypes. Unknown pairs are not verified negatives. A future evaluation must remove target/reverse/alias equivalents and prevent indication text or related relations from trivially revealing the answer; disease-family or source/temporal holdouts are preferable to a naive row split.

A lexical scan for `colorectal`, whole-word `colon`, or `rectal` yields 85 disease candidates, including colorectal cancer, colon adenocarcinoma, colorectal adenocarcinoma, rectal cancer, benign lesions, hereditary categories, and possible irrelevant name matches. This is an inclusive audit list, **not a curated equivalence mapping**. Across these candidates there are 23 indication pairs and 12 drugs; see `primekg_colorectal_disease_candidates.csv` and `primekg_crc_indications.csv` for exact records. GI.CIN, GI.HM-indel, GI.HM-SNV and GI.GS must not be merged into those broad disease nodes.

The key modeling gap is **cross-node-type transfer**: a model trained only on Drug→Disease indications is not automatically trained to score Drug→Subtype edges. New subtype nodes have no therapeutic labels; an unseen relation decoder or unsupported mapping to Disease cannot be assumed to work. A shared disease/subtype representation, independently justified subtype-of links, or an explicit transfer objective could be considered later, but each requires validation against independent subtype-sensitive evidence. No technical compatibility statistic solves this gap.

## Minimal integration proposal (not implemented)

| Proposed relation | Meaning and constraints |
| --- | --- |
| Sample → has observed gene flag → Gene | Repaired source observations; preserve sample provenance and binary semantics |
| Sample → has subtype → CRCSubtype | Existing assignments; retain missing assignments |
| CRCSubtype → has observed gene → Gene | New aggregation relation with numerator, denominator, frequency and source; do not call it causal or treats |
| Gene → participates in → SourcePathway | ECoGraph membership, retaining collection and unverified identities |
| Gene ↔ PrimeKG gene/protein | One entity only after approved unambiguous ID resolution; retain all original aliases and provenance |
| Drug → targets → Gene/protein | Existing PrimeKG `drug_protein` + `target`, distinct from enzyme/carrier/transporter |
| Drug → indicated for → Disease | Existing PrimeKG `indication`, distinct from contraindications/off-label records |

Gene matches can bridge subtype observations into a drug-containing molecular graph. Preserve CRCSubtype as a new semantic node type and introduce a new aggregate observation relation; do not silently reuse Disease or a therapeutic relation. ECoGraph named pathways should remain namespaced nodes until an audited stable-ID crosswalk exists, otherwise duplicate/equivalent pathway entities will be created. PrimeKG genes sharing the same NCBI ID also need a deliberate canonicalization policy. With these changes, heterogeneous message passing is technically possible; meaningful learning and response prediction remain hypotheses.

## Recommended next steps

1. Regenerate and validate sample–gene edges from the local matrix, explicitly excluding patient metadata; supply pathway vertices. Review the exact binary mutation preprocessing and patient-level sampling policy. The current source export itself is **not currently merge-ready**.
2. Obtain an authoritative symbol→stable gene ID crosswalk, review unmatched and ambiguous symbols, and obtain a versioned MSigDB/Reactome pathway crosswalk. Preserve rather than guess unverified equivalents.
3. Define subtype molecular representations using per-sample frequencies, variant consequences if obtainable, and burden-controlled associations. Validate apparent subtype differences in an independent cohort; six GI.HM-SNV samples limit confidence.
4. Define the prediction target and transfer design before training. Establish independent subtype-sensitive therapeutic evidence for evaluation. Compare popularity, direct-target proximity and pathway/network baselines; a GNN must add measured value beyond them.
5. Only after these checks, consider graph integration and separate GNN experiments. Keep an eventual investigative agent downstream of validated graph/model tools and mark predictions as hypotheses.

## Reproduction and scope

Run `python3 analysis/run_compatibility.py` from the project root. `analysis/README.md` documents dependencies and optional suppression of supplementary recovery. The script never modifies datasets or existing preprocessing outputs. PrimeKG is processed in chunks with only relevant drug edges retained; no second full edge list or merged graph is saved. All original mandatory inputs are SHA-256 checked before and after execution. `reports/compatibility/README.md` indexes artifacts. Statistical screening results are supplementary, and no drugs are recommended as treatments.
