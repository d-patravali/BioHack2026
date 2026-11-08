"""Summarize a completed phenotype GNN run without training or modifying sources."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def table(frame):
    def render(value):
        if isinstance(value,(float,np.floating)):return f'{value:.4f}'
        return str(value)
    return '\n'.join(['| '+' | '.join(frame.columns)+' |','| '+' | '.join(['---']*len(frame.columns))+' |',*[ '| '+' | '.join(render(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)]])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path)
    args=parser.parse_args();run=args.run
    manifest=json.loads((run/'manifest.json').read_text())
    if not manifest.get('input_hashes_unchanged'): raise ValueError('Run has no completed preservation audit')
    metrics=pd.read_csv(run/'metrics.csv')
    risk=pd.read_csv(run/'selective_prediction.csv')
    splits=pd.read_csv(run/'splits.csv')
    intact=metrics[metrics.scenario.eq('intact')]
    overview=intact.groupby(['method','condition'])[['mrr','gene_macro_mrr','hits1','hits10','hits100']].mean().reset_index()
    sensitivity=metrics[metrics.method.eq('relational_gnn')].groupby(['condition','scenario'])[['mrr','hits10','hits100']].mean().reset_index()
    selective=risk.groupby(['condition','desired_validation_coverage'])[['test_coverage','accepted_gene_macro_mrr','all_gene_macro_mrr']].mean().reset_index()
    pair_metrics = pd.read_csv(run/'per_pair_metrics.csv.gz')
    coverage = pd.read_csv(Path(__file__).resolve().parents[1]/'reports/inference/gene_coverage.csv')
    pair_metrics = pair_metrics.merge(coverage[['gene_index','ppi_neighbors','pathways']],on='gene_index',validate='many_to_one')
    strata = []
    for dimension in ['ppi_neighbors','pathways']:
        staged = pair_metrics[pair_metrics.scenario.eq('intact')].copy()
        staged['stratum'] = np.where(staged[dimension].gt(0),'present','absent')
        grouped = staged.groupby(['seed','method','condition','stratum']).agg(mrr=('mrr','mean'),hits10=('hits10','mean'),pairs=('mrr','size')).reset_index()
        grouped['dimension'] = dimension
        strata.append(grouped)
    pd.concat(strata).to_csv(run/'context_strata.csv',index=False)
    cold=[]
    for seed,rows in splits.groupby('seed'):
        train=rows[rows.split.eq('train')];test=rows[rows.split.eq('test')];val=rows[rows.split.eq('validation')]
        assert not set(train.gene_index)&set(test.gene_index)
        assert not set(train.gene_index)&set(val.gene_index)
        assert not set(val.gene_index)&set(test.gene_index)
        cold.append({'seed':seed,'train_genes':train.gene_index.nunique(),'validation_genes':val.gene_index.nunique(),'test_genes':test.gene_index.nunique(),'train_pairs':len(train),'validation_pairs':len(val),'test_pairs':len(test),'test_pairs_without_training_phenotype':int((~test.phenotype_index.isin(train.phenotype_index)).sum())})
    neural=overview[overview.method.eq('relational_gnn')]
    pathway=overview.loc[overview.method.eq('shared_pathway_transfer'),'mrr'].iloc[0]
    pop=overview.loc[overview.method.eq('phenotype_popularity'),'mrr'].iloc[0]
    best=neural.mrr.max()
    conclusion=('A fitted GNN exceeds both pathway transfer and popularity in mean MRR on these masks.' if best>max(pathway,pop)
                else 'The tested GNN does not exceed the strongest heuristic in mean MRR. This run does not establish a learned-model advantage.')
    checkpoints=metrics[(metrics.method.eq('relational_gnn')) & metrics.scenario.eq('intact')][['seed','condition','best_epoch','validation_mrr']]
    text=f'''# Gene–phenotype prediction under missing evidence

## Outcome

**{conclusion}** This is an executed exploratory benchmark, not a new biological discovery or an independently validated novelty claim. The CRC app, notebooks, historical results, and datasets were not modified by this experiment.

## Target and protocol

Predict explicitly recorded `phenotype_protein` associations for genes whose **entire direct phenotype annotation set is withheld**. There are {manifest['canonical_relation_pairs']['phenotype_protein']:,} canonical positive pairs; candidates are all {manifest['node_counts']['effect/phenotype']:,} phenotype nodes. Disease paths are never promoted to positive labels.

The previous gene-holdout test masks (seeds {manifest['settings']['seeds']}) are reused and checked exactly by original canonical node IDs. A further 20% of annotated genes is reserved for validation, leaving approximately 60/20/20 training/validation/test genes. Earlier results used 80/20, so their numerical values are not a directly matched comparison; all heuristics here are recomputed with the same reduced supervision as the GNN. These test masks were already inspected historically: this is exploratory evaluation, not an untouched confirmatory test.

{table(pd.DataFrame(cold))}

All query genes remain in the molecular graph: this is **transductive annotation-cold evaluation**, not unseen-node generalization or a test on naturally unannotated genes. Three masks reuse overlapping genes and one fixed initialization seed (42); means are descriptive, not confidence intervals.

## Model and leakage controls

A small PyTorch encoder uses typed learned ID embeddings (width {manifest['settings']['width']}), two relation-specific mean-message layers with residual state, and a gene–phenotype dot-product decoder with phenotype bias. Training uses balanced positive and sampled-unlabeled logistic terms (five unlabeled draws per positive). Unlabeled draws exclude training positives; they are optimization contrasts, not verified biological negatives. No held-out labels enter sampling or message passing.

- **annotations_only:** training gene–phenotype edges and their reverse messages, without PPI/pathways/disease context; a learned-ID/annotation control.
- **molecular:** adds PPI and gene–pathway relations in both directions where applicable.
- **with_disease:** also adds disease–gene and positive disease–phenotype relations. Shared underlying evidence can make this channel semantically leaky even after target removal.

Exact `(type, source, id)` aliases are canonicalized before splitting. Reciprocal rows are verified and deduplicated. Held-out target edges and reverse copies are absent from the training adjacency. Checkpoints are selected using validation MRR only, evaluated every {manifest['settings']['evaluate_every']} epochs up to {manifest['settings']['epochs']}. No phenotype hierarchy/name features, patient data, drug labels, or CRC subtype nodes are used.

{table(checkpoints)}

## Matched intact-evidence results

Mean over masks; per-mask and per-pair values are saved separately. MRR averages reciprocal target rank; Hits@K is the expected fraction recovered within K ranks. Other recorded positives for each query gene are filtered **only during evaluation**. Exact score ties use the existing expected-rank implementation, avoiding arbitrary ordering. Missing associations remain unknown, so these measure annotation recovery, not biological specificity.

{table(overview)}

The strongest comparison is against pathway transfer and fixed PPI/pathway fusion, not merely a zero-score lookup. Popularity ignores the query gene. The learned annotation-only control tests whether context adds value beyond the encoder's learned IDs and label prior.

The matched pathway baseline has mean MRR **{pathway:.4f}**, compared with popularity **{pop:.4f}** and the best mean GNN condition **{best:.4f}**. The molecular GNN and the disease-enabled GNN do not establish useful incremental graph learning in this run. Context-availability strata (PPI/pathway present versus absent) are also saved in `context_strata.csv`; they are descriptive and not matched-degree causal comparisons.

## Evidence-removal experiments

Each selected checkpoint is evaluated after removing an entire relation channel, including its reverse, **at inference without retraining or checkpoint reselection**. This tests distribution-shift sensitivity, not the causal truth or independent value of a biological mechanism. Unavailable channels produce an unchanged control. `with_disease` versus `molecular` provides a separately retrained comparison with/without disease channels, although other architectural details/parameters differ.

{table(sensitivity)}

Heuristics receive the analogous channel removals (`metrics.csv`). Removing training annotations from the GNN leaves learned parameters intact: it is not erasing every trace of training knowledge. It is a structural sensitivity test.

## Exploratory abstention

Confidence is the unfiltered top-1 minus top-2 score margin. Thresholds are chosen from **validation genes only**, at desired retained fractions 100%, 75%, 50%, and 25%, and transferred unchanged to test genes. All threshold ties are retained, so actual test coverage can differ substantially. No score margin is called a calibrated probability. We compare accepted-gene macro MRR with all-gene macro MRR, and include identical diagnostics for the heuristics.

{table(selective)}

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
python3 analysis/phenotype_report.py {run.as_posix()}
python3 -m unittest discover -s tests/research -p 'test_phenotype_gnn.py' -v
```

Dependencies: existing research requirements plus PyTorch (CPU). No PyTorch Geometric, external API, downloads, or agent credentials required. The runner creates a new timestamped folder; an existing explicit `--output` is rejected. The report generator writes only this run's `findings.md` and summary tables.

Artifacts: `metrics.csv`, `per_pair_metrics.csv.gz`, `splits.csv`, `training_curves.csv`, `selective_prediction.csv`, `top_hypotheses.csv`, `context_strata.csv`, `model_*.pt`, and `manifest.json`. Hypotheses are bounded top-five examples; arbitrary within-tie order is not meaningful. The manifest records settings, versions, input/code hashes and preservation checks. Runtime: {manifest['elapsed_seconds']:.1f} seconds. Protected input hashes unchanged: **{manifest['input_hashes_unchanged']}**.
'''
    (run/'findings.md').write_text(text)
    overview.to_csv(run/'summary.csv',index=False)
    sensitivity.to_csv(run/'robustness_summary.csv',index=False)
    print(conclusion)
    print(run/'findings.md')

if __name__=='__main__':main()
