"""Audit a completed CMS drug-scoring run and build an unverified review queue."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score
from crc_report import markdown,tie_top_mask


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('run',type=Path);args=parser.parse_args();run=args.run
    manifest=json.loads((run/'manifest.json').read_text());assert manifest['input_hashes_unchanged']
    metrics=pd.read_csv(run/'evaluation_metrics.csv')
    summary=metrics.groupby('condition').agg(runs=('test_ap','size'),mean_ap=('test_ap','mean'),min_ap=('test_ap','min'),max_ap=('test_ap','max'),mean_disease_macro_ap=('test_disease_macro_ap','mean')).reset_index()
    diagnostic=pd.read_csv(run/'subtype_pairwise_diagnostics.csv')
    sensitivity=diagnostic[diagnostic.profile.isin(['signed','magnitude','permuted_genes'])].groupby(['condition','profile'])[['tie_inclusive_jaccard','score_spearman','score_rms_difference']].mean().reset_index()
    variance=pd.read_csv(run/'score_variance_decomposition.csv')
    vsummary=variance[variance.profile.isin(['signed','magnitude'])].groupby(['condition','profile']).drug_main_effect_variance_pct.mean().reset_index()
    ties=pd.read_csv(run/'tie_audit.csv')
    tied=ties[ties.profile.isin(['signed','magnitude'])].groupby(['condition','profile']).agg(min_top_tie=('top_score_tie_count','min'),max_top_tie=('top_score_tie_count','max'),min_inclusive_top20=('inclusive_top20_count','min'),max_inclusive_top20=('inclusive_top20_count','max')).reset_index()
    candidates=pd.read_csv(run/'candidate_hypotheses.csv.gz',dtype={'recorded_target_ncbi_ids':str,'drug_id':str})
    consensus=candidates.groupby(['condition','profile','subtype','drug_id','drug_name','recorded_target_ncbi_ids'],dropna=False).agg(seeds_in_top20=('seed','nunique'),best_tied_rank=('tied_rank_min','min'),mean_logit=('raw_logit','mean'),std_logit=('raw_logit','std'),mean_within_drug_contrast=('within_drug_subtype_contrast','mean'),training_indication_count=('training_indication_count','first')).reset_index()
    overlay=pd.read_csv(run/'cms_expression_edges.csv',dtype={'y_id':str})
    signature_sets={cms:set(overlay.loc[overlay.x_id.str.endswith(':'+cms),'y_id']) for cms in ['CMS1','CMS2','CMS3','CMS4']}
    consensus['direct_signature_target_overlap']=consensus.apply(lambda r:len(set(str(r.recorded_target_ncbi_ids).split('|'))&signature_sets[r.subtype]),axis=1)
    consensus['model_seeds_total']=len(manifest['initialization_seeds'])
    consensus['status']='exploratory scoring consistency, not biological validation'
    consensus.to_csv(run/'candidate_consensus.csv',index=False)
    # The pathway variant is prespecified, not selected after test performance.
    # Keep signed and magnitude transfer assumptions separate; never label novelty.
    queue=consensus[(consensus.condition.eq('ppi_pathway')) & consensus.profile.isin(['signed','magnitude'])].copy()
    queue=queue.sort_values(['profile','subtype','seeds_in_top20','best_tied_rank','drug_id'],ascending=[True,True,False,True,True],kind='stable')
    queue['search_query']=queue.apply(lambda r:f'"{r.drug_name}" colorectal cancer "{r.subtype}"',axis=1)
    for field in ['supporting_evidence','contradictory_evidence','prior_crc_use','experimental_model','dose_and_mechanism','citation','reviewer_conclusion']:
        queue[field]='NOT REVIEWED'
    queue['clinical_status']='No therapeutic efficacy or novelty established'
    queue.to_csv(run/'literature_review_queue.csv',index=False)
    evidence=run/'test_predictions.npz';arr=np.load(evidence)
    disease_groups=[np.flatnonzero(arr['disease_positions']==d) for d in np.unique(arr['disease_positions'])]
    bootstrap=[];rng=np.random.default_rng(2026)
    for repeat in range(500):
        take=np.concatenate([disease_groups[i] for i in rng.choice(len(disease_groups),len(disease_groups),replace=True)])
        y=arr['labels'][take];base=average_precision_score(y,arr['drug_popularity'][take])
        for condition in ['ppi','ppi_pathway','drug_only']:
            for seed in manifest['initialization_seeds']:
                ap=average_precision_score(y,arr[f'{condition}_{seed}'][take])
                bootstrap.append({'repeat':repeat,'condition':condition,'seed':seed,'delta_ap_vs_popularity':ap-base})
    ci=pd.DataFrame(bootstrap).groupby(['condition','seed']).delta_ap_vs_popularity.quantile([.025,.5,.975]).unstack().reset_index()
    ci.columns=['condition','seed','lower','median','upper'];ci.to_csv(run/'paired_disease_bootstrap.csv',index=False)
    baseline=float(summary.loc[summary.condition.eq('drug_popularity'),'mean_ap'].iloc[0])
    modelmeans=summary[~summary.condition.eq('drug_popularity')]
    outcome='Popularity exceeds every fitted GNN in this run.' if metrics[~metrics.condition.eq('drug_popularity')].test_ap.max()<baseline else 'Consult the full matched metrics: the learned results are not uniformly below popularity.'
    schedules=[]
    for profile in ['signed','magnitude']:
        subset=candidates[candidates.condition.eq('ppi_pathway')&candidates.profile.eq(profile)]
        for (seed,subtype),rows in subset.groupby(['seed','subtype']):
            top=rows[rows.tied_rank_min.eq(1)]
            schedules.append({'profile':profile,'seed':seed,'subtype':subtype,'leading_drugs':len(top),'distinct_leading_target_signatures':top.recorded_target_ncbi_ids.nunique(),'leading_target_signatures':' ; '.join(sorted(top.recorded_target_ncbi_ids.dropna().unique()))})
    leading=pd.DataFrame(schedules);leading.to_csv(run/'leading_target_audit.csv',index=False)
    leading_targets='; '.join(sorted(leading.leading_target_signatures.unique()))
    overlap_range=queue.direct_signature_target_overlap.agg(['min','max'])
    signature_sets=run/'profile_summary.csv'
    data_summary=manifest['graph']
    report=f'''# CRCSC–PrimeKG GNN drug-hypothesis prototype

## Verdict

**Technical integration and exploratory scoring executed. {outcome} Meaningful subtype-specific drug prioritization remains unvalidated.** Do not present these scores as treatment recommendations, efficacy probabilities, or literature-novel discoveries.

The full PrimeKG CSV is referenced read-only. A small typed CMS overlay is exported and compiled into the runtime molecular graph. This avoids copying the 8.1M-row graph while preserving the original broader KG for future evidence retrieval. The encoder uses a selected gene/PPI/pathway/drug/disease view; phenotype and other PrimeKG relations do not enter this model. The working agent remains pinned to its earlier reviewed CRC backend; it has not silently switched to CMS results.

## Integration and semantics

- Four new `crc_cms_subtype` nodes (CMS1–CMS4) bridge through **400 discovery-selected expression edges** to existing canonical NCBI gene nodes.
- `higher_expression_signature` and `lower_expression_signature` remain distinct relations. Both signed weight and magnitude, adjusted discovery effect, FDR and source are saved. We do not reinterpret them as mutations, causal dependencies, or therapeutic indications.
- Exact `(type, source, id)` aliases are canonicalized before training. The base-node audit preserves original indices and canonical representatives. The compact graph has {data_summary['genes']:,} genes, {data_summary['drugs_all_roles']:,} drug nodes in the matched evaluation universe and {data_summary['drugs_with_direct_targets']:,} direct-target-supported hypothesis candidates.
- PPI: {data_summary['ppi_directed_pairs']:,} directed messages; gene–pathway membership: {data_summary['gene_pathway_pairs']:,}; disease–gene associations: {data_summary['disease_gene_pairs']:,}. The full graph is not relabeled or overwritten.
- Only `drug_protein` display role `target` enters drug pooling. Enzymes, carriers and transporters are not treated as targets. Reciprocal target copies and exact-ID aliases are deduplicated.

{markdown(pd.read_csv(signature_sets))}

The source representation is `{manifest['source_run']}`. Discovery-only profile weights are frozen; validation replication statistics were not used to select drugs or fit model parameters. CMS labels remain distinct from ECoGraph GI genomic subtypes.

## Runtime message-passing view

```mermaid
flowchart LR
    CMS[CMS subtype] -->|higher/lower expression context| Gene[PrimeKG gene]
    Gene -->|PPI| OtherGene[PrimeKG gene]
    Gene -->|membership| Pathway[Reactome pathway]
    Drug[PrimeKG drug] -->|recorded target| Gene
    Disease[PrimeKG disease] -->|associated gene| Gene
    Drug -.->|indication supervision only| Disease
    CMS -.->|unvalidated score| Drug
```

## Model and the key transfer assumption

We reuse the existing width48 PPI/target/disease pooling architecture and 40-epoch-style training recipe, adding one relation-specific pathway propagation variant. **ppi** uses PPI; **ppi_pathway** adds gene→pathway→gene mean messages; **drug_only** removes molecular disease context from the decoder while retaining drug target/PPI representations. Pathway pooling includes a gene's own pathway membership. These are compact message-passing models, not an end-to-end encoder of all PrimeKG relation families.

Known disease indications train the decoder; no CMS drug labels exist. Therapeutic relations never enter molecular adjacency. At inference, signed CMS weights pool gene embeddings into the shared disease projection. **This transfers a decoder trained on unsigned disease associations to signed expression-derived subtype vectors. That domain-transfer assumption is not validated.** A negative expression weight is not an instruction to inhibit or activate a gene. Signed scoring is exploratory, not a learned pharmacological-direction model.

Prespecified input variants are signed, magnitude-only, higher-expression-only, and lower-expression-only. CMS4 has no lower-expression genes in the selected top100; its empty lower-only profile is a zero-input diagnostic, and no candidate hypotheses are exported from that empty row. Weighted profile changes occur at inference and cannot improve disease-indication test AP in this one-way pooling architecture.

## Matched held-out disease evaluation

All conditions share a seed42 disease split, indication positives, fixed validation/test unknowns, and epoch-specific training unknown pairs. Checkpoint selection uses validation AP only. Split details: `{manifest['split']}`. Sampled unknowns exclude recorded indications, contraindications and off-label relationships; they are not confirmed biological negatives. Disease-family/source overlap remains possible. Candidate drugs remain identical across conditions, including drugs without direct targets for matched evaluation; subtype hypotheses are restricted to direct-target-supported drugs.

{markdown(summary)}

Average precision is sampled association recovery, not accuracy or clinical benefit. Three initialization seeds on a single split are not independent cohort replication. Every selected checkpoint should also be checked against the epoch budget in the raw metrics; boundary checkpoints do not establish convergence. Canonical gene indexing differs from the older experiment, so model scores need not numerically reproduce that run; the paired methods within this run have the same inputs/splits. `paired_disease_bootstrap.csv` resamples test diseases to describe AP differences against popularity, conditional on these labels, sampling scheme and fitted models. It does not establish source-independent validation.

## CMS sensitivity and exact ties

**Observed leading group:** the pathway model's signed and magnitude profiles have {leading.leading_drugs.min()}–{leading.leading_drugs.max()} exactly tied leading drugs across all CMS/seed combinations. Their recorded target signatures are `{leading_targets}` (NCBI2908 is NR3C1). Direct overlap with the selected subtype signatures ranges from {int(overlap_range['min'])} to {int(overlap_range['max'])} target genes. Zero direct overlap does not exclude indirect graph connections, but these results do not establish subtype-specific therapeutic support. Cross-seed persistence can reflect a stable architectural shortcut rather than reliable biology.

These are inference diagnostics, not therapeutic performance metrics. Top20 sets retain the **entire exact cutoff tie**. Jaccard measures overlap of those inclusive sets; raw logit separation is not efficacy separation.

{markdown(sensitivity)}

{markdown(tied)}

Drug-level main-effect percentage of total score variance:

{markdown(vsummary)}

The remaining variance includes subtype offsets and drug-by-subtype interaction; this decomposition is descriptive, not causal. Identical recorded target neighborhoods necessarily produce identical embeddings and scores in this architecture, which lacks compound-specific input features. More selective CMS inputs cannot distinguish those drugs by themselves. See `leading_target_audit.csv` for the actual leading target signatures rather than interpreting arbitrary within-tie ordering.

Zero and identical profiles produce identical subtype score vectors. Swapping profiles swaps prediction vectors (asserted). Gene permutation is a null input: altered scores or rankings do not show biological correctness. The drug-only model must ignore subtype profiles (asserted). `profile_perturbations.csv` and full `scores_*.npz` retain these diagnostics. Within-drug subtype contrasts subtract each drug's four-subtype mean for inspection only; they are not a validated alternative predictor and were not used to select the primary hypotheses.

## Candidate files and downstream literature review

`candidate_hypotheses.csv.gz` contains tie-inclusive top20 drug scores for every actual nonempty profile, fitted condition, seed, and subtype. It includes DrugBank ID, name, raw logit, minimum tied rank, full tie-group size, training-indication popularity, and recorded target signatures. Cross-seed consistency is in `candidate_consensus.csv`; it is not a probability of correctness. Scores from different fitted models are not calibrated or directly comparable in magnitude. Direct signature-target overlap counts are saved with the consensus/queue; zero overlap is absence of a direct recorded connection, not absence of an indirect path or a proof of biological irrelevance.

`literature_review_queue.csv` uses the **prespecified pathway variant** under signed and magnitude profiles, without selecting a model because its test result looks better. It preserves all inclusive candidate groups and leaves evidence fields **NOT REVIEWED**. Sort order prioritizes recurrence across seeds for investigation, not efficacy. **No live literature review has been performed.** Searches should assess:

1. Whether the compound/target association and mechanism are correct, including pharmacological direction and off-target effects.
2. Evidence specifically concerning CRC/CMS versus broad cancer results, and cell-line/animal/patient context.
3. Supporting, null, and contradictory results, dose/exposure, toxicity, and potential tissue-composition confounding.
4. Whether an association is already known; missing graph links or search results do not establish novelty.

Literature may contextualize or weaken a model hypothesis; it does not independently verify clinical efficacy. A genuinely subtype-sensitive ranking needs independent response/dependency labels and a transfer design evaluated against simple baselines. Genomic/microenvironment context, signed expression, and drug action should not be conflated.

## Reproduction and artifacts

```sh
python3 analysis/crcsc_drug_gnn.py
python3 analysis/crcsc_drug_report.py {run.as_posix()}
python3 -m unittest discover -s tests/research -p 'test_crcsc_drug_gnn.py' -v
```

Research requirements plus PyTorch, scikit-learn and threadpoolctl; no API key or download required. New output directories only. The report generator reads completed results and writes report/audit/queue files, without retraining or changing the original experiment. Saved checkpoints require the same audited graph ordering for reconstruction.

Files include `cms_nodes.csv`, `cms_expression_edges.csv`, `base_node_alias_audit.csv.gz`, `profile_gene_audit.csv`, `profile_summary.csv`, `gene_positions.csv`, `disease_split.csv`, `evaluation_pairs.npz`, `evaluation_metrics.csv`, `test_predictions.npz`, `training_curves.csv`, `model_*.pt`, `scores_*.npz`, `candidate_hypotheses.csv.gz`, `candidate_consensus.csv`, `literature_review_queue.csv`, `leading_target_audit.csv`, `tie_audit.csv`, `profile_perturbations.csv`, `subtype_pairwise_diagnostics.csv`, `score_variance_decomposition.csv`, `paired_disease_bootstrap.csv`, and `manifest.json`.

Protected source/representation/notebook/agent/UI hashes unchanged: **{manifest['input_hashes_unchanged']}**. This is a prototype integration and hypothesis audit; the output has not been promoted to the agent backend.
'''
    (run/'findings.md').write_text(report)
    summary.to_csv(run/'evaluation_summary.csv',index=False)
    print(summary.to_string(index=False));print('Report:',run/'findings.md');print('Review queue:',run/'literature_review_queue.csv')

if __name__=='__main__':main()
