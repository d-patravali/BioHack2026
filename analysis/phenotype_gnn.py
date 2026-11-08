"""Exploratory relation-aware gene–phenotype completion under missing evidence.

Original data/results stay read-only. A new output directory is mandatory.
"""
import argparse
import hashlib
import json
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from torch import nn
from torch.nn import functional as F

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from explore_inference import expected_metrics

RELATIONS = {
    'phenotype_protein': ('gene/protein', 'effect/phenotype'),
    'protein_protein': ('gene/protein', 'gene/protein'),
    'pathway_protein': ('gene/protein', 'pathway'),
    'disease_protein': ('gene/protein', 'disease'),
    'disease_phenotype_positive': ('disease', 'effect/phenotype'),
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def load_graph():
    nodes = pd.read_csv(ROOT / 'data/processed/primekg/nodes.csv', keep_default_na=False,
                        dtype={'id': str}).sort_values('index')
    canonical = nodes.groupby(['type', 'source', 'id'])['index'].transform('min')
    mapping = dict(zip(nodes['index'], canonical))
    representatives = nodes[nodes['index'].eq(canonical)].copy()
    types = ['gene/protein', 'effect/phenotype', 'pathway', 'disease']
    ids = {kind: representatives.loc[representatives.type.eq(kind), 'index'].to_numpy() for kind in types}
    local = {kind: {int(v): i for i, v in enumerate(values)} for kind, values in ids.items()}
    edges = {r: set() for r in RELATIONS}
    reverse = {r: set() for r in RELATIONS}
    usecols = ['x_index', 'y_index', 'x_type', 'y_type', 'relation']
    for chunk in pd.read_csv(ROOT / 'data/raw/primekg/kg.csv', usecols=usecols, chunksize=250000):
        for rel, (left, right) in RELATIONS.items():
            part = chunk[chunk.relation.eq(rel)]
            for x, y, xt, yt in part[['x_index', 'y_index', 'x_type', 'y_type']].itertuples(index=False, name=None):
                x, y = int(mapping[x]), int(mapping[y])
                if xt == left and yt == right:
                    edges[rel].add((local[left][x], local[right][y]))
                if xt == right and yt == left:
                    reverse[rel].add((local[left][y], local[right][x]))
    for rel in edges:
        assert edges[rel] == reverse[rel], f'Reciprocal mismatch: {rel}'
    # Canonical aliases can create self loops even if the raw file has none.
    edges['protein_protein'] = {(a,b) for a,b in edges['protein_protein'] if a != b}
    return ids, edges, representatives


def split_labels(labels, seed):
    annotated = np.array(sorted({g for g,p in labels}))
    rng = np.random.default_rng(seed)
    test = set(rng.choice(annotated, round(.2 * len(annotated)), replace=False).tolist())
    remaining = np.array(sorted(set(annotated) - test))
    validation = set(np.random.default_rng(seed + 1000).choice(
        remaining, round(.2 * len(annotated)), replace=False).tolist())
    groups = {'test': test, 'validation': validation, 'train': set(remaining) - validation}
    result = {k: sorted((g,p) for g,p in labels if g in genes) for k,genes in groups.items()}
    assert not groups['train'] & groups['validation'] and not groups['train'] & groups['test']
    assert not groups['validation'] & groups['test']
    assert set().union(*(set(v) for v in result.values())) == labels
    return result


def row_normalized(pairs, shape):
    pairs = sorted(pairs)
    if not pairs:
        return sp.csr_matrix(shape, dtype=np.float32)
    a,b = np.array(pairs).T
    matrix = sp.csr_matrix((np.ones(len(a), dtype=np.float32), (a,b)), shape=shape)
    return sp.diags(1 / np.maximum(1, np.asarray(matrix.sum(1)).ravel())) @ matrix


def as_torch(matrix):
    matrix = matrix.tocoo()
    return torch.sparse_coo_tensor(np.array([matrix.row, matrix.col]), matrix.data.astype(np.float32), matrix.shape).coalesce()


def graph_adjacencies(ids, edges, train, condition):
    n = {k: len(v) for k,v in ids.items()}
    specs = {'annotations': ('gene/protein', 'effect/phenotype', set(train))}
    if condition != 'annotations_only':
        specs['ppi'] = ('gene/protein', 'gene/protein', edges['protein_protein'])
        specs['pathways'] = ('gene/protein', 'pathway', edges['pathway_protein'])
    if condition == 'with_disease':
        specs['disease_genes'] = ('gene/protein', 'disease', edges['disease_protein'])
        specs['disease_phenotypes'] = ('disease', 'effect/phenotype', edges['disease_phenotype_positive'])
    result = {}
    for name,(left,right,pairs) in specs.items():
        result[name] = (left, right, as_torch(row_normalized(pairs, (n[left], n[right]))))
        if left != right:
            result[name + '_reverse'] = (right, left, as_torch(row_normalized({(b,a) for a,b in pairs}, (n[right],n[left]))))
    return result


class RelationalGNN(nn.Module):
    """Two relation-specific mean-message layers, residual state and dot decoder.

    Typed learned ID embeddings make this a transductive annotation-cold model,
    not an unseen-node or feature-based inductive model.
    """
    def __init__(self, sizes, adjacency, width=24):
        super().__init__()
        self.kinds = list(sizes)
        self.keys = {k: 't'+str(i) for i,k in enumerate(self.kinds)}
        self.embeddings = nn.ParameterDict({self.keys[k]: nn.Parameter(torch.randn(v,width)*.1) for k,v in sizes.items()})
        self.adjacency = adjacency
        self.layers = nn.ModuleList([nn.ModuleDict({name: nn.Linear(width,width,bias=False) for name in adjacency}) for _ in range(2)])
        self.self_layers = nn.ModuleList([nn.ModuleDict({self.keys[k]: nn.Linear(width,width) for k in sizes}) for _ in range(2)])
        self.bias = nn.Parameter(torch.zeros(sizes['effect/phenotype']))
        self.scale = width ** .5

    def encode(self, removed=()):
        h = {k: self.embeddings[self.keys[k]] for k in self.kinds}
        for layer, own in zip(self.layers, self.self_layers):
            result = {k: own[self.keys[k]](h[k]) for k in self.kinds}
            degree = Counter()
            for name,(left,right,adj) in self.adjacency.items():
                if name.removesuffix('_reverse') in removed:
                    continue
                result[left] = result[left] + layer[name](torch.sparse.mm(adj,h[right]))
                degree[left] += 1
            h = {k: F.relu(result[k] / max(1, degree[k]+1)) + .2*h[k] for k in self.kinds}
        return h

    def score_pairs(self, h, pairs):
        return (h['gene/protein'][pairs[:,0]]*h['effect/phenotype'][pairs[:,1]]).sum(1)/self.scale + self.bias[pairs[:,1]]

    def score_all(self, h, genes):
        return (h['gene/protein'][genes] @ h['effect/phenotype'].T / self.scale + self.bias).cpu().numpy()


def sample_unknown(train, nphenotypes, rng, per_positive=5):
    """Unlabeled contrasts, not verified negatives; do not consult held-out labels."""
    forbidden = set(train)
    result = []
    for g,p in train:
        for _ in range(per_positive):
            candidate = int(rng.integers(nphenotypes))
            while (g,candidate) in forbidden:
                candidate = int(rng.integers(nphenotypes))
            result.append((g,candidate))
    assert not set(result) & forbidden
    return np.asarray(result)


def evaluate(scores, queries, truth, all_labels, nphenotypes):
    harmonic = np.r_[0., np.cumsum(1 / np.arange(1,nphenotypes+1))]
    full = {}
    for g,p in all_labels:
        full.setdefault(g,set()).add(p)
    rows = []
    by_gene = {g: scores[i] for i,g in enumerate(queries)}
    for g,p in truth:
        values = expected_metrics(by_gene[g], p, full[g], harmonic)
        rows.append({'gene': g, 'phenotype': p, **dict(zip(['mrr','hits1','hits10','hits100'],values))})
    return pd.DataFrame(rows)


def heuristic_scores(ids, edges, train, queries, removed=()):
    ng,nphe = len(ids['gene/protein']),len(ids['effect/phenotype'])
    # Unnormalized labels, exactly the same training supervision as the GNN.
    y = sp.csr_matrix((np.ones(len(train)),np.array(train).T),shape=(ng,nphe))
    ppi = row_normalized(edges['protein_protein'],(ng,ng))
    # Existing baseline counts PPI votes, not mean votes.
    ppi = sp.diags(np.diff(ppi.indptr).astype(np.float32)) @ ppi
    memberships = edges['pathway_protein']
    a,b=np.array(sorted(memberships)).T
    membership=sp.csr_matrix((np.ones(len(a)),(a,b)),shape=(ng,len(ids['pathway'])))
    count=np.asarray(membership.sum(0)).ravel()
    path = membership @ sp.diags(1/np.maximum(1,count-1)) @ (membership.T @ y)
    # Query genes are annotation-cold; their own y rows are zero.
    ppi_scores = (ppi @ y)[queries].toarray()
    path_scores = path[queries].toarray()
    if 'ppi' in removed: ppi_scores *= 0
    if 'pathways' in removed: path_scores *= 0
    popularity = np.broadcast_to(np.asarray(y.sum(0)),(len(queries),nphe))
    gd=row_normalized(edges['disease_protein'],(ng,len(ids['disease'])))
    gd=sp.diags(np.diff(gd.indptr).astype(np.float32)) @ gd
    dp=row_normalized(edges['disease_phenotype_positive'],(len(ids['disease']),nphe))
    disease=(gd @ dp)[queries].toarray()
    if 'disease_genes' in removed or 'disease_phenotypes' in removed: disease *= 0
    return {'phenotype_popularity': popularity,'ppi_transfer':ppi_scores,
            'shared_pathway_transfer':path_scores,'disease_paths':disease,
            'equal_ppi_pathway': .5*(ppi_scores/np.maximum(1,ppi_scores.max(1,keepdims=True)) + path_scores/np.maximum(1,path_scores.max(1,keepdims=True)))}


def selective_rows(seed, condition, scores, validation_scores, queries, frame):
    """Validation-only margin cutoffs; coverage includes all threshold ties."""
    ordered = np.sort(validation_scores, axis=1)
    valmargin = ordered[:, -1] - ordered[:, -2]
    ordered = np.sort(scores, axis=1)
    margin = ordered[:, -1] - ordered[:, -2]
    gene_mrr = frame.groupby('gene').mrr.mean()
    result = []
    for retained in [1., .75, .5, .25]:
        threshold = float(np.quantile(valmargin, 1-retained)) if retained < 1 else -float('inf')
        accepted = [g for g,m in zip(queries, margin) if m >= threshold]
        result.append({'seed': seed, 'condition': condition,
            'desired_validation_coverage': retained, 'margin_threshold': threshold,
            'test_coverage': len(accepted)/len(queries), 'accepted_genes': len(accepted),
            'accepted_gene_macro_mrr': float(gene_mrr.reindex(accepted).mean()) if accepted else None,
            'all_gene_macro_mrr': float(gene_mrr.mean())})
    return result


def fit(ids, edges, split, condition, args, out, seed):
    torch.manual_seed(42)
    model=RelationalGNN({k:len(v) for k,v in ids.items()},graph_adjacencies(ids,edges,split['train'],condition),args.width)
    optimizer=torch.optim.AdamW(model.parameters(),lr=.01,weight_decay=.001)
    positive=torch.tensor(split['train'])
    rng=np.random.default_rng(seed+2000)
    validation_genes=sorted({g for g,p in split['validation']})
    labels=set().union(*(set(v) for v in split.values()))
    best=-1; history=[]; saved=None
    for epoch in range(1,args.epochs+1):
        model.train();optimizer.zero_grad();h=model.encode()
        unknown=torch.tensor(sample_unknown(split['train'],len(ids['effect/phenotype']),rng))
        pos=model.score_pairs(h,positive);neg=model.score_pairs(h,unknown)
        loss=.5*(F.softplus(-pos).mean()+F.softplus(neg).mean())
        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
        row={'split_seed':seed,'condition':condition,'epoch':epoch,'loss':float(loss.detach())}
        if epoch % args.evaluate_every==0 or epoch==args.epochs:
            model.eval()
            with torch.no_grad():
                scores=model.score_all(model.encode(),validation_genes)
            metric=evaluate(scores,validation_genes,split['validation'],labels,len(ids['effect/phenotype'])).mrr.mean()
            row['validation_mrr']=metric
            if metric>best:
                best=metric;saved={k:v.detach().clone() for k,v in model.state_dict().items()};best_epoch=epoch
        history.append(row)
    model.load_state_dict(saved);model.eval()
    torch.save({'state_dict':saved,'condition':condition,'width':args.width,'best_epoch':best_epoch,'split_seed':seed},out/f'model_{condition}_{seed}.pt')
    return model,history,best_epoch,best


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--epochs',type=int,default=60)
    parser.add_argument('--seeds',type=int,nargs='+',default=[17,29,43])
    parser.add_argument('--width',type=int,default=24)
    parser.add_argument('--evaluate-every',type=int,default=10)
    args=parser.parse_args()
    if min(args.epochs,args.width,args.evaluate_every)<1: parser.error('Positive training settings required')
    out=args.output or ROOT/'reports/phenotype_gnn'/datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(2)
    start=time.monotonic()
    protected=[ROOT/'data/raw/primekg/kg.csv',ROOT/'data/processed/primekg/nodes.csv',ROOT/'reports/inference/test_splits.csv']
    hashes={str(p.relative_to(ROOT)):digest(p) for p in protected}
    print(f'Loading graph; output={out}',flush=True)
    ids,edges,nodes=load_graph()
    labels=edges['phenotype_protein']
    assert len(labels)==3330
    old=pd.read_csv(ROOT/'reports/inference/test_splits.csv')
    summary={'node_counts':{k:len(v) for k,v in ids.items()},'canonical_relation_pairs':{k:len(v) for k,v in edges.items()},'settings':vars(args)|{'output':str(out)},'initialization_seed':42,'annotation_holdout':'60/20/20 by annotated gene; existing test masks reused; exploratory not untouched tests'}
    (out/'manifest.json').write_text(json.dumps(summary,indent=2))
    metrics=[];details=[];curves=[];splitrows=[];riskrows=[];prediction_rows=[]
    for seed in args.seeds:
        split=split_labels(labels,seed)
        # Check exact biological IDs, not arbitrary local-index ordering.
        previous=old[(old['mode']=='gene_annotation_holdout') & old.seed.eq(seed)]
        actual={(int(ids['gene/protein'][g]),int(ids['effect/phenotype'][p])) for g,p in split['test']}
        assert actual==set(zip(previous.gene_index,previous.phenotype_index))
        for name,pairs in split.items():
            splitrows.extend({'seed':seed,'split':name,'gene_index':int(ids['gene/protein'][g]),'phenotype_index':int(ids['effect/phenotype'][p])} for g,p in pairs)
        queries=sorted({g for g,p in split['test']})
        def record(method,scenario,scores,condition='heuristic'):
            frame=evaluate(scores,queries,split['test'],labels,len(ids['effect/phenotype']))
            frame['seed']=seed;frame['method']=method;frame['scenario']=scenario;frame['condition']=condition
            # Audit external IDs; local IDs retained to allow joins within the run.
            frame['gene_index']=frame.gene.map(lambda g:int(ids['gene/protein'][g]))
            frame['phenotype_index']=frame.phenotype.map(lambda p:int(ids['effect/phenotype'][p]))
            details.append(frame)
            row={'seed':seed,'method':method,'scenario':scenario,'condition':condition,'test_pairs':len(frame),'test_genes':len(queries)}
            row.update(frame[['mrr','hits1','hits10','hits100']].mean().to_dict())
            row['gene_macro_mrr']=frame.groupby('gene').mrr.mean().mean();metrics.append(row)
            return frame
        valgenes = sorted({g for g,p in split['validation']})
        valheuristics = heuristic_scores(ids, edges, split['train'], valgenes)
        for scenario,removed in [('intact',()),('no_ppi',('ppi',)),('no_pathways',('pathways',)),('no_disease',('disease_genes','disease_phenotypes'))]:
            for method,scores in heuristic_scores(ids,edges,split['train'],queries,removed).items():
                frame = record(method,scenario,scores)
                if scenario == 'intact':
                    riskrows.extend(selective_rows(seed,method,scores,valheuristics[method],queries,frame))
        del valheuristics
        for condition in ['annotations_only','molecular','with_disease']:
            print(f'Split {seed}: fitting {condition}',flush=True)
            model,history,best_epoch,best=fit(ids,edges,split,condition,args,out,seed)
            curves.extend(history)
            valgenes=sorted({g for g,p in split['validation']})
            with torch.no_grad():
                valscores=model.score_all(model.encode(),valgenes)
            for scenario,removed in [('intact',()),('no_ppi',('ppi',)),('no_pathways',('pathways',)),('no_disease',('disease_genes','disease_phenotypes')),('no_training_annotations',('annotations',))]:
                with torch.no_grad(): scores=model.score_all(model.encode(removed),queries)
                frame=record('relational_gnn',scenario,scores,condition)
                metrics[-1].update(best_epoch=best_epoch,validation_mrr=best)
                if scenario=='intact':
                    riskrows.extend(selective_rows(seed,condition,scores,valscores,queries,frame))
                    # Save a bounded set of examples; unknown high scores are hypotheses.
                    for i,g in enumerate(queries):
                        for p in np.argsort(-scores[i],kind='stable')[:5]:
                            prediction_rows.append({'seed':seed,'condition':condition,'gene_index':int(ids['gene/protein'][g]),'phenotype_index':int(ids['effect/phenotype'][p]),'score':float(scores[i,p]),'recorded_test_pair':(g,int(p)) in set(split['test'])})
            pd.DataFrame(metrics).to_csv(out/'metrics.csv',index=False)
            pd.DataFrame(curves).to_csv(out/'training_curves.csv',index=False)
            print(f'  best epoch {best_epoch}, validation MRR {best:.4f}; results saved',flush=True)
    pd.concat(details).to_csv(out/'per_pair_metrics.csv.gz',index=False,compression='gzip')
    pd.DataFrame(splitrows).to_csv(out/'splits.csv',index=False)
    pd.DataFrame(riskrows).to_csv(out/'selective_prediction.csv',index=False)
    pd.DataFrame(prediction_rows).to_csv(out/'top_hypotheses.csv',index=False)
    summary['elapsed_seconds']=time.monotonic()-start
    summary['input_sha256']=hashes
    summary['input_hashes_unchanged']=all(digest(ROOT/p)==v for p,v in hashes.items())
    summary['code_sha256']={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),ROOT/'scripts/explore_inference.py']}
    summary['versions']={'torch':torch.__version__,'numpy':np.__version__,'pandas':pd.__version__}
    assert summary['input_hashes_unchanged']
    (out/'manifest.json').write_text(json.dumps(summary,indent=2))
    print(pd.DataFrame(metrics).query("scenario == 'intact'").groupby(['method','condition'])[['mrr','hits10','hits100']].mean().to_string(),flush=True)
    print(f'Completed: {out}',flush=True)


if __name__=='__main__':
    main()
