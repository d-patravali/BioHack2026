"""CRCSC overlay + compact PrimeKG GNN hypothesis experiment; never efficacy.

Keeps the full PrimeKG as a read-only referenced base, exports a small typed
CMS overlay, and trains a molecular message-passing view with matched controls.
"""
import argparse
import json
import sys
from copy import deepcopy
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from crc_experiments import IndicationGNN, weighted_adjacency, make_split, pair_metrics, fingerprint
from crc_report import tie_top_mask, markdown

CMS=['CMS1','CMS2','CMS3','CMS4']
LABELS=['indication','off-label use','contraindication']


def normalize_signed(weights):
    weights=np.asarray(weights,dtype=np.float32)
    if weights.ndim!=2 or not np.isfinite(weights).all():raise ValueError("Finite two-dimensional profile matrix required")
    denominator=np.abs(weights).sum(1,keepdims=True)
    return np.divide(weights,denominator,out=np.zeros_like(weights),where=denominator>0)


def prepare(signatures, out):
    original=pd.read_csv(ROOT/'data/processed/primekg/nodes.csv',dtype={'id':str},keep_default_na=False)
    canonical=original.groupby(['type','source','id'])['index'].transform('min')
    aliases=dict(zip(original['index'],canonical))
    nodes=original[original['index'].eq(canonical)].sort_values('index').set_index('index')
    genes=pd.Index(nodes.index[nodes.type.eq('gene/protein')])
    chunks=[]
    relations=['protein_protein','pathway_protein','drug_protein','disease_protein',*LABELS]
    for chunk in pd.read_csv(ROOT/'data/raw/primekg/kg.csv',usecols=['x_index','y_index','relation','display_relation'],chunksize=250000):
        part=chunk[chunk.relation.isin(relations)].copy()
        part['x_index']=part.x_index.map(aliases);part['y_index']=part.y_index.map(aliases)
        chunks.append(part)
    edges=pd.concat(chunks,ignore_index=True)
    def pairs(relation,left,right,role=None):
        rows=edges[edges.relation.isin(relation)]
        if role:rows=rows[rows.display_relation.eq(role)]
        forward=rows.x_index.map(nodes.type).eq(left)
        result=pd.DataFrame({'left':rows.x_index.where(forward,rows.y_index),'right':rows.y_index.where(forward,rows.x_index)}).drop_duplicates()
        assert result.left.map(nodes.type).eq(left).all() and result.right.map(nodes.type).eq(right).all()
        return result
    ppi=pairs(['protein_protein'],'gene/protein','gene/protein')
    ppi=ppi[ppi.left.ne(ppi.right)]
    all_roles=pairs(['drug_protein'],'drug','gene/protein')
    direct=pairs(['drug_protein'],'drug','gene/protein','target')
    gd=pairs(['disease_protein'],'disease','gene/protein')
    gp=pairs(['pathway_protein'],'gene/protein','pathway')
    known=pairs(['indication'],'disease','drug')
    excluded=pairs(LABELS,'disease','drug')
    drugs=pd.Index(sorted(all_roles.left.unique()))
    subtype_indices=np.arange(int(original['index'].max())+1,int(original['index'].max())+5)
    diseases=pd.Index(sorted(set(gd.left)|set(subtype_indices)))
    paths=pd.Index(sorted(gp.right.unique()))
    gpos=pd.Series(np.arange(len(genes)),index=genes)
    dpos=pd.Series(np.arange(len(diseases)),index=diseases)
    rpos=pd.Series(np.arange(len(drugs)),index=drugs)
    ppos=pd.Series(np.arange(len(paths)),index=paths)
    adjacency=lambda left,right,nleft,nright:weighted_adjacency(left,right,np.ones(len(left)),nleft,nright)
    adj_ppi=adjacency(gpos.loc[ppi.left],gpos.loc[ppi.right],len(genes),len(genes))
    adj_drug=adjacency(rpos.loc[direct.left],gpos.loc[direct.right],len(drugs),len(genes))
    adj_gpath=adjacency(gpos.loc[gp.left],ppos.loc[gp.right],len(genes),len(paths))
    adj_pathg=adjacency(ppos.loc[gp.right],gpos.loc[gp.left],len(paths),len(genes))
    signatures=signatures.copy()
    signatures['ncbi_id']=signatures.ncbi_id.astype(str)
    id_to_node=nodes[nodes.type.eq('gene/protein')].reset_index().set_index('id')['index']
    assert id_to_node.index.is_unique
    signatures['gene_node_index']=signatures.ncbi_id.map(id_to_node)
    assert signatures.gene_node_index.notna().all()
    assert not signatures.duplicated(['subtype','ncbi_id']).any()
    signed=np.zeros((4,len(genes)),dtype=np.float32)
    overlay=[];newnodes=[]
    for i,cms in enumerate(CMS):
        identifier='crcsc:GSE39582:'+cms
        newnodes.append({'node_index':int(subtype_indices[i]),'node_id':identifier,'node_type':'crc_cms_subtype','node_name':cms,'source':'CRCSC GSE39582 discovery expression signature','definition':'expression-based CMS; not equivalent to ECoGraph GI genomic subtype'})
        for row in signatures[signatures.subtype.eq(cms)].itertuples(index=False):
            index=int(row.gene_node_index);value=float(row.profile_weight)
            signed[i,gpos.loc[index]]=value
            overlay.append({'x_index':int(subtype_indices[i]),'x_id':identifier,'x_type':'crc_cms_subtype','relation':'higher_expression_signature' if value>0 else 'lower_expression_signature','y_index':index,'y_id':row.ncbi_id,'y_namespace':'NCBI','y_type':'gene/protein','signed_weight':value,'magnitude':abs(value),'discovery_adjusted_effect':row.discovery_adjusted_effect,'discovery_q':row.discovery_q,'source':'discovery only; validation evidence not used for profile selection','interpretation':'relative expression context, not causal dependency or drug action direction'})
    signed=normalize_signed(signed)
    pd.DataFrame([{'subtype':cms,'genes':int((signed[i]!=0).sum()),'higher_expression_genes':int((signed[i]>0).sum()),'lower_expression_genes':int((signed[i]<0).sum()),'l1_mass':float(np.abs(signed[i]).sum())} for i,cms in enumerate(CMS)]).to_csv(out/'profile_summary.csv',index=False)
    profiles={'signed':signed,'magnitude':normalize_signed(np.abs(signed)),'up_only':normalize_signed(np.maximum(signed,0)),'down_only':normalize_signed(np.maximum(-signed,0))}
    assert np.allclose(np.abs(signed).sum(1),1)
    base=adjacency(dpos.loc[gd.left],gpos.loc[gd.right],len(diseases),len(genes)).coalesce()
    si,sj=np.nonzero(signed)
    # This computational reduction is an explicit transfer assumption: the
    # disease decoder was supervised on unsigned disease-gene associations.
    ii=np.vstack([dpos.loc[subtype_indices].to_numpy()[si],sj])
    adj_disease=torch.sparse_coo_tensor(torch.cat([base.indices(),torch.tensor(ii)],dim=1),torch.cat([base.values(),torch.tensor(signed[si,sj])]),base.shape).coalesce()
    positives=known[known.left.isin(diseases)&known.right.isin(drugs)].copy()
    target_mask=drugs.isin(direct.left)
    pd.DataFrame(newnodes).to_csv(out/'cms_nodes.csv',index=False)
    pd.DataFrame(overlay).to_csv(out/'cms_expression_edges.csv',index=False)
    original[['index','type','source','id']].assign(canonical_index=canonical).to_csv(out/'base_node_alias_audit.csv.gz',index=False,compression='gzip')
    signatures.to_csv(out/'profile_gene_audit.csv',index=False)
    summary={'genes':len(genes),'diseases_including_cms':len(diseases),'drugs_all_roles':len(drugs),'drugs_with_direct_targets':int(target_mask.sum()),'ppi_directed_pairs':len(ppi),'gene_pathway_pairs':len(gp),'disease_gene_pairs':len(gd),'eligible_indications':len(positives),'cms_nodes':4,'cms_signature_edges':len(overlay),'training_uses_expression_sign':False,'cms_therapeutic_labels':0}
    # Existing helpers consume these keys; no ECoGraph inputs are loaded.
    return {'nodes':nodes,'genes':genes,'drugs':drugs,'diseases':diseases,'positives':positives,'excluded':excluded,'target_mask':target_mask,'adj_ppi':adj_ppi,'adj_drug':adj_drug,'adj_disease':adj_disease,'adj_gpath':adj_gpath,'adj_pathg':adj_pathg,'profiles':profiles,'direct':direct,'summary':summary}


class CMSDrugGNN(IndicationGNN):
    def __init__(self,data,use_pathways=False,use_context=True):
        super().__init__(len(data['genes']),data['adj_ppi'],data['adj_drug'],data['adj_disease'],use_context=use_context)
        self.use_pathways=use_pathways
        self.adj_gpath,self.adj_pathg=data['adj_gpath'],data['adj_pathg']
        self.gene_path=nn.Linear(48,48,bias=False) if use_pathways else None

    def encode(self):
        initial=self.gene.weight
        messages=self.gene_self(initial)+self.gene_ppi(torch.sparse.mm(self.adj_ppi,initial))
        if self.use_pathways:
            path=torch.sparse.mm(self.adj_pathg,initial)
            messages=messages+self.gene_path(torch.sparse.mm(self.adj_gpath,path))
        gene_vectors=F.relu(messages)
        drugs=self.drug_projection(torch.sparse.mm(self.adj_drug,gene_vectors))
        pooled=torch.sparse.mm(self.adj_disease,gene_vectors)
        if not self.use_context:pooled=torch.zeros_like(pooled)
        return self.disease_projection(pooled),drugs,gene_vectors


def fit(data,split,condition,seed,epochs,out):
    positive,negative,train_unknown,_=split
    torch.manual_seed(seed)
    model=CMSDrugGNN(data,use_pathways=condition=='ppi_pathway',use_context=condition!='drug_only')
    optimizer=torch.optim.AdamW(model.parameters(),lr=.003,weight_decay=.001)
    best=-1;state=None;curve=[]
    for epoch in range(epochs):
        model.train();optimizer.zero_grad()
        pairs=np.concatenate([positive['train'],train_unknown[epoch]])
        y=torch.tensor(np.r_[np.ones(len(positive['train'])),np.zeros(len(train_unknown[epoch]))],dtype=torch.float32)
        loss=F.binary_cross_entropy_with_logits(model.score(model.encode(),pairs),y)
        loss.backward();optimizer.step();model.eval()
        with torch.no_grad():
            vp=np.concatenate([positive['validation'],negative['validation']])
            vy=np.r_[np.ones(len(positive['validation'])),np.zeros(len(negative['validation']))]
            ap=average_precision_score(vy,model.score(model.encode(),vp).numpy())
        curve.append({'condition':condition,'seed':seed,'epoch':epoch+1,'loss':float(loss.detach()),'validation_ap':ap})
        if ap>best:best=ap;best_epoch=epoch+1;state=deepcopy(model.state_dict())
    model.load_state_dict(state);model.eval()
    metrics,scores=pair_metrics(model,positive['test'],negative['test'],data['target_mask'])
    row={'condition':condition,'seed':seed,'best_epoch':best_epoch,'validation_ap':best,'test_ap':metrics['ap'],'test_disease_macro_ap':metrics['disease_macro_ap'],'test_target_supported_ap':metrics['target_supported_ap']}
    torch.save({'state_dict':state,'condition':condition,'seed':seed,'width':48},out/f'model_{condition}_{seed}.pt')
    print(json.dumps(row),flush=True)
    return model,row,curve,scores


def score_profiles(model,data):
    profiles=dict(data['profiles']);primary=profiles['signed']
    profiles['zero']=np.zeros_like(primary)
    profiles['identical']=np.repeat(primary.mean(0,keepdims=True),4,axis=0)
    profiles['swapped']=primary[[1,0,3,2]]
    profiles['permuted_genes']=primary[:,np.random.default_rng(2026).permutation(primary.shape[1])]
    result={}
    with torch.no_grad():
        _,drug,gene=model.encode();positions=np.flatnonzero(data['target_mask'])
        for name,weights in profiles.items():
            pooled=torch.tensor(weights)@gene
            if not model.use_context:pooled=torch.zeros_like(pooled)
            subtype=model.disease_projection(pooled)
            result[name]=np.stack([model.score_vectors(s.expand(len(positions),-1),drug[positions]).numpy() for s in subtype])
    np.testing.assert_allclose(result['swapped'],result['signed'][[1,0,3,2]],atol=1e-5)
    np.testing.assert_allclose(result['zero'],np.repeat(result['zero'][:1],4,axis=0),atol=1e-5)
    np.testing.assert_allclose(result['identical'],np.repeat(result['identical'][:1],4,axis=0),atol=1e-5)
    if not model.use_context:np.testing.assert_allclose(result['signed'],np.repeat(result['signed'][:1],4,axis=0),atol=1e-5)
    return result,positions


def audit_rankings(arrays,positions,data,condition,seed):
    drugs=data['drugs'][positions];metadata=data['nodes']
    direct=data['direct'];targetsets=direct.groupby('left').right.agg(lambda x:sorted(set(x)))
    signatures={int(d):'|'.join(sorted(set(metadata.loc[gs,'id']))) for d,gs in targetsets.items()}
    candidates=[];diagnostics=[];perturb=[];ties=[];variance=[]
    primary=arrays['signed']
    for profile,scores in arrays.items():
        total=np.square(scores-scores.mean()).sum();drug_effect=np.square(scores.mean(0)-scores.mean()).sum()*4
        variance.append({'condition':condition,'seed':seed,'profile':profile,'drug_main_effect_variance_pct':100*drug_effect/total if total else None})
        masks=[tie_top_mask(s,20) for s in scores]
        for i,j in combinations(range(4),2):
            intersection=int((masks[i]&masks[j]).sum());union=int((masks[i]|masks[j]).sum())
            diagnostics.append({'condition':condition,'seed':seed,'profile':profile,'subtype_a':CMS[i],'subtype_b':CMS[j],'tie_inclusive_intersection':intersection,'tie_inclusive_union':union,'tie_inclusive_jaccard':intersection/union,'score_spearman':float(spearmanr(scores[i],scores[j]).statistic),'score_rms_difference':float(np.sqrt(np.mean((scores[i]-scores[j])**2)))})
        for i,cms in enumerate(CMS):
            leading=scores[i]==scores[i].max()
            ties.append({'condition':condition,'seed':seed,'profile':profile,'subtype':cms,'top_score_tie_count':int(leading.sum()),'inclusive_top20_count':int(masks[i].sum())})
            perturb.append({'condition':condition,'seed':seed,'profile':profile,'subtype':cms,'rms_change_vs_signed':float(np.sqrt(np.mean((scores[i]-primary[i])**2))),'top20_jaccard_vs_signed':float((masks[i]&tie_top_mask(primary[i],20)).sum()/(masks[i]|tie_top_mask(primary[i],20)).sum())})
            if profile not in data['profiles'] or np.abs(data['profiles'][profile][i]).sum()==0:continue
            for p in np.flatnonzero(masks[i]):
                node=int(drugs[p]);value=scores[i,p]
                candidates.append({'condition':condition,'seed':seed,'profile':profile,'subtype':cms,'drug_node_index':node,'drug_id':metadata.loc[node,'id'],'drug_name':metadata.loc[node,'name'],'tied_rank_min':int((scores[i]>value).sum())+1,'tie_group_size':int((scores[i]==value).sum()),'raw_logit':float(value),'training_indication_count':int(data['training_indication_counts'][data['drugs'].get_loc(node)]),'within_drug_subtype_contrast':float(value-scores[:,p].mean()),'recorded_target_ncbi_ids':signatures[node],'status':'unvalidated subtype-drug association hypothesis; no efficacy, novelty, or drug-action direction established'})
    return candidates,diagnostics,perturb,ties,variance


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run',type=Path,default=ROOT/'reports/crcsc/20261108_172638')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--epochs',type=int,default=40)
    parser.add_argument('--seeds',type=int,nargs='+',default=[42,43,44])
    args=parser.parse_args()
    if args.epochs<1:parser.error('Positive epoch count required')
    out=args.output or ROOT/'reports/crcsc_drug_gnn'/datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    out.mkdir(parents=True,exist_ok=False)
    protected=[ROOT/'data/raw/primekg/kg.csv',ROOT/'data/processed/primekg/nodes.csv',*args.source_run.glob('*')]
    protected += [p for folder in ['notebooks','agent','ui'] for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in str(p)]
    hashes={str(p.relative_to(ROOT)):fingerprint(p) for p in protected if p.is_file()}
    print('Preparing canonical graph and CMS overlay:',out,flush=True)
    signatures=pd.read_csv(args.source_run/'selected_gene_signatures.csv',dtype={'ncbi_id':str})
    source=json.loads((args.source_run/'manifest.json').read_text());assert source['input_hashes_unchanged']
    data=prepare(signatures,out)
    print(json.dumps(data['summary']),flush=True)
    torch.set_num_threads(2);torch.use_deterministic_algorithms(True)
    split=make_split(data,out,args.epochs,split_seed=42)
    positive,negative,_,ids=split
    pairs=np.concatenate([positive['test'],negative['test']]);y=np.r_[np.ones(len(positive['test'])),np.zeros(len(negative['test']))]
    counts=np.bincount(positive['train'][:,1],minlength=len(data['drugs']))
    data['training_indication_counts']=counts
    supported=data['target_mask'][pairs[:,1]]
    baseline={'condition':'drug_popularity','seed':42,'best_epoch':0,'validation_ap':None,'test_ap':average_precision_score(y,counts[pairs[:,1]]),'test_disease_macro_ap':np.mean([average_precision_score(y[pairs[:,0]==d],counts[pairs[pairs[:,0]==d,1]]) for d in np.unique(pairs[:,0])]),'test_target_supported_ap':average_precision_score(y[supported],counts[pairs[supported,1]])}
    metrics=[baseline];curves=[];candidates=[];diagnostics=[];perturb=[];ties=[];variance=[]
    evals={'labels':y,'disease_positions':pairs[:,0],'drug_positions':pairs[:,1],'drug_popularity':counts[pairs[:,1]]}
    for seed in args.seeds:
        for condition in ['ppi','ppi_pathway','drug_only']:
            model,row,curve,scores=fit(data,split,condition,seed,args.epochs,out)
            metrics.append(row);curves.extend(curve);evals[f'{condition}_{seed}']=scores
            arrays,positions=score_profiles(model,data)
            np.savez_compressed(out/f'scores_{condition}_{seed}.npz',drug_node_indices=data['drugs'][positions].to_numpy(),**arrays)
            a,b,c,d,e=audit_rankings(arrays,positions,data,condition,seed)
            candidates.extend(a);diagnostics.extend(b);perturb.extend(c);ties.extend(d);variance.extend(e)
            pd.DataFrame(metrics).to_csv(out/'evaluation_metrics.csv',index=False)
            pd.DataFrame(curves).to_csv(out/'training_curves.csv',index=False)
    pd.DataFrame(candidates).to_csv(out/'candidate_hypotheses.csv.gz',index=False,compression='gzip')
    pd.DataFrame(diagnostics).to_csv(out/'subtype_pairwise_diagnostics.csv',index=False)
    pd.DataFrame(perturb).to_csv(out/'profile_perturbations.csv',index=False)
    pd.DataFrame(ties).to_csv(out/'tie_audit.csv',index=False)
    pd.DataFrame(variance).to_csv(out/'score_variance_decomposition.csv',index=False)
    np.savez_compressed(out/'test_predictions.npz',**evals)
    pd.DataFrame({'gene_position':np.arange(len(data['genes'])),'gene_node_index':data['genes'],'ncbi_id':data['nodes'].loc[data['genes'],'id'].to_numpy()}).to_csv(out/'gene_positions.csv',index=False)
    manifest={'source_run':str(args.source_run.relative_to(ROOT)),'graph':data['summary'],'epochs':args.epochs,'initialization_seeds':args.seeds,'split_seed':42,'split':{name:{'diseases':len(ds),'positives':len(positive[name]),'unknowns':len(negative[name]) if name!='train' else len(positive[name])} for name,ds in ids.items()},'candidate_universe':'all direct-target-supported PrimeKG drugs; full cutoff ties retained','training_target':'known disease indications versus sampled unknown pairs','cms_training_labels':0,'phenotype_edges_used':False,'signed_transfer_assumption':'CMS signed expression pooling is transferred through a decoder trained on unsigned disease-gene associations; action direction and transfer validity are unverified','input_sha256':hashes,'input_hashes_unchanged':all(fingerprint(ROOT/p)==v for p,v in hashes.items()),'code_sha256':{str(p.relative_to(ROOT)):fingerprint(p) for p in [Path(__file__),ROOT/'analysis/crc_experiments.py',ROOT/'analysis/crc_report.py']},'torch_version':torch.__version__,'numpy_version':np.__version__,'pandas_version':pd.__version__}
    assert manifest['input_hashes_unchanged']
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(pd.DataFrame(metrics).groupby('condition').test_ap.mean().to_string(),flush=True)
    print('Complete:',out,flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
