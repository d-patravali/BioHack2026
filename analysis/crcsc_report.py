"""Produce a concise interpretation of a completed CRCSC experiment; no retraining."""
import argparse
import json
from pathlib import Path
import pandas as pd


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('run',type=Path);args=parser.parse_args();run=args.run
    manifest=json.loads((run/'manifest.json').read_text())
    assert manifest['input_hashes_unchanged']
    metrics=pd.read_csv(run/'classification_metrics.csv').set_index('model')
    replication=pd.read_csv(run/'profile_replication.csv')
    coverage=pd.read_csv(run/'signature_drug_connectivity.csv')
    ci=pd.read_csv(run/'paired_bootstrap_accuracy_delta.csv').set_index('model').loc['pathway_logistic']
    genes=replication[replication.level.eq('gene')]
    pathway=replication[replication.level.eq('pathway')]
    shuffled=metrics.loc[[f'shuffled_pathway_{s}' for s in [17,29,43]],'accuracy']
    text=f'''# CRCSC assessment: a stronger molecular starting point, not validated drug inference

**Branch:** `CRCSC-branch`. **Completed run:** `{run.name}`.

## What was established

- All 566 expression samples join to CMS labels and clinical records. The final labels provide 519 classified cases (CMS1:91, CMS2:232, CMS3:69, CMS4:127); 47 unclassified cases are excluded from supervised analysis.
- The supplied split yields 409 classified discovery and 110 classified validation cases. Features, scaling, signature selection, and classifier fitting use discovery only.
- {manifest['single_id_probes']:,} probes have unambiguous single-NCBI-ID mappings, yielding {manifest['unique_annotated_genes']:,} genes. **{manifest['primekg_mapped_genes']:,} map to PrimeKG ({manifest['primekg_mapping_pct']:.2f}%)**. This provides stable-ID bridges without fuzzy gene-name matching.
- Each CMS receives 100 discovery-selected genes with signed adjusted expression weights. Direction agreement in validation is {100*genes.selected_sign_agreement.min():.0f}–{100*genes.selected_sign_agreement.max():.0f}%; **{genes.selected_validation_q05_same_direction.min()}–{genes.selected_validation_q05_same_direction.max()} of 100** also satisfy validation BH q<=.05 in the same direction. Batch/age/sex/stage/location adjustment is descriptive and does not establish causality.
- {manifest['pathways_min5_genes']:,} PrimeKG Reactome pathways have at least five measured genes. Across all pathway effects, discovery–validation Spearman correlations range **{pathway.all_effect_spearman.min():.3f}–{pathway.all_effect_spearman.max():.3f}**. Selected pathway direction agreement is {100*pathway.selected_sign_agreement.min():.0f}–{100*pathway.selected_sign_agreement.max():.0f}%.

## Simple classifiers and graph controls

| Diagnostic | Validation accuracy | Balanced accuracy |
| --- | ---: | ---: |
| Majority | {100*metrics.loc['majority','accuracy']:.1f}% | {100*metrics.loc['majority','balanced_accuracy']:.1f}% |
| Batch only | {100*metrics.loc['batch_only','accuracy']:.1f}% | {100*metrics.loc['batch_only','balanced_accuracy']:.1f}% |
| Expression logistic | {100*metrics.loc['expression_logistic','accuracy']:.1f}% | {100*metrics.loc['expression_logistic','balanced_accuracy']:.1f}% |
| Signed 100-gene signature scores | {100*metrics.loc['signed_signature_scores','accuracy']:.1f}% | {100*metrics.loc['signed_signature_scores','balanced_accuracy']:.1f}% |
| Pathway logistic | {100*metrics.loc['pathway_logistic','accuracy']:.1f}% | {100*metrics.loc['pathway_logistic','balanced_accuracy']:.1f}% |

Shuffled-pathway accuracy ranges **{100*shuffled.min():.1f}–{100*shuffled.max():.1f}%**. The pathway-versus-expression accuracy difference is {100*(metrics.loc['pathway_logistic','accuracy']-metrics.loc['expression_logistic','accuracy']):.2f} percentage points; its paired within-class bootstrap 95% interval is **[{100*ci.lower:.2f}, {100*ci.upper:.2f}] percentage points**, including zero. This panel does not establish a reliable advantage from curated pathway membership. The bootstrap captures uncertainty within this validation set, not upstream label uncertainty or external-cohort generalization.

CMS labels were derived from expression, and GSE39582 contributed to the original CMS work. High reproduction accuracy is expected to some degree and must not be pitched as new subtype discovery, independent biological validation, or a demonstrated need for a GNN.

## Connection to the intended subtype–drug goal

The goal remains a model that uses subtype molecular context in a larger graph to score **subtype–drug hypotheses**. This experiment tests the reliability and compatibility of that context before reusing it in the model. It does not train a GNN or drug-response predictor.

The 100-gene signatures have the following recorded direct-target connectivity:

| Subtype | Signature genes with direct drug targets | Unique connected drugs |
| --- | ---: | ---: |
'''
    for row in coverage.itertuples(index=False):text+=f'| {row.subtype} | {row.selected_genes_with_targets}/100 | {row.unique_connected_drugs} |\n'
    text+='''
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
'''
    (run/'assessment.md').write_text(text)
    print(run/'assessment.md')

if __name__=='__main__':main()
