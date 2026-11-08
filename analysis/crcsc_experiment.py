"""Read-only CRCSC molecular-profile reproducibility experiment.

Uses supplied discovery/validation groups, not drug labels. Creates a new run.
No graph merge, GNN fitting, or change to the existing application.
"""
import argparse
import gzip
import hashlib
import io
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import sparse, stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, confusion_matrix
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.pipeline import make_pipeline
import sklearn
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
CMS=['CMS1','CMS2','CMS3','CMS4']
FINAL='CMS_final_network_plus_RFclassifier_in_nonconsensus_samples'
ANNOTATION_URL='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(2**20),b''):h.update(block)
    return h.hexdigest()


def bh(p):
    p=np.asarray(p,dtype=float);order=np.argsort(p);ranked=p[order]
    q=np.minimum.accumulate((ranked*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=np.empty_like(q);out[order]=np.clip(q,0,1)
    return out


def parse_gene_id(value):
    """Accept exactly one distinct numeric Entrez identifier; never pick an alias."""
    ids=set(re.findall(r'\d+',str(value)))
    residual=re.sub(r'\d+|[/;|\s]', '',str(value))
    if str(value).strip() in {'','NA','---','nan'}:return '', 'unmapped'
    if residual:return '', 'invalid_identifier'
    if len(ids)!=1:return '', 'ambiguous'
    return next(iter(ids)), 'single_ncbi_id'


def read_annotation(path):
    header=[];body=[];inside=False
    with gzip.open(path,'rt') as stream:
        for line in stream:
            if line.startswith('!platform_table_begin'):inside=True;continue
            if line.startswith('!platform_table_end'):break
            if inside:body.append(line)
            else:header.append(line.rstrip())
    frame=pd.read_csv(io.StringIO(''.join(body)),sep='\t',dtype=str,keep_default_na=False)
    assert {'ID','Gene ID','Gene symbol'}<=set(frame)
    assert frame.ID.is_unique
    parsed=[parse_gene_id(x) for x in frame['Gene ID']]
    frame['ncbi_id']=[p[0] for p in parsed];frame['mapping_status']=[p[1] for p in parsed]
    return frame[['ID','Gene ID','Gene symbol','ncbi_id','mapping_status']],header


def sample_id(value):
    found=re.match(r'^(GSM\d+)(?:_|$)',value)
    if not found:raise ValueError(f'Unexpected expression sample identifier: {value}')
    return found.group(1)


def load_expression(path, annotation):
    with path.open() as stream:columns=stream.readline().rstrip('\r\n').split('\t')
    samples=[sample_id(c) for c in columns]
    assert len(samples)==len(set(samples))
    # This file's header omits a name for the first/probe-ID field.
    frame=pd.read_csv(path,sep='\t',header=None,skiprows=1,names=['probe_id']+samples,index_col=0,dtype={s:np.float32 for s in samples})
    assert frame.index.is_unique and np.isfinite(frame.to_numpy()).all()
    audit=pd.DataFrame({'probe_id':frame.index}).merge(annotation,left_on='probe_id',right_on='ID',how='left',validate='1:1')
    audit['mapping_status']=audit.mapping_status.fillna('missing_annotation')
    eligible=audit[audit.mapping_status.eq('single_ncbi_id')]
    selected=frame.loc[eligible.probe_id].copy();selected.index=eligible.ncbi_id.to_numpy()
    # Fixed within-sample aggregation, without selecting a probe by subtype/validation.
    genes=selected.groupby(level=0,sort=True).median()
    assert genes.index.is_unique
    return genes,audit,columns


def standardized(train, test):
    mean=train.mean(0);scale=train.std(0,ddof=1)
    valid=np.isfinite(scale)&(scale>1e-6)
    return (train[:,valid]-mean[valid])/scale[valid],(test[:,valid]-mean[valid])/scale[valid],valid


def covariates(frame):
    # Nuisance design is fitted separately within each cohort for adjusted
    # descriptive contrasts, never used to transform discovery from validation.
    age=pd.to_numeric(frame.age,errors='coerce');age=age.fillna(age.median())
    age=(age-age.mean())/max(float(age.std()),1)
    categories=frame[['batch','gender','stage','tumorLocation']].astype(str)
    dummy=pd.get_dummies(categories,drop_first=True,dtype=float)
    return np.column_stack([np.ones(len(frame)),age.to_numpy(),dummy.to_numpy()])


def adjusted_contrast(values, group, design):
    """OLS subtype-vs-rest contrast controlling clinical/batch covariates."""
    design=np.asarray(design,dtype=float);group=np.asarray(group,dtype=float)
    projected=np.linalg.pinv(design)
    resid_group=group-design@(projected@group)
    resid_y=values-design@(projected@values)
    denominator=resid_group@resid_group
    df=len(group)-np.linalg.matrix_rank(np.column_stack([design,group]))
    if denominator<1e-8 or df<1:raise ValueError('Subtype contrast not identifiable after covariate adjustment')
    beta=resid_group@resid_y/denominator
    residual=resid_y-resid_group[:,None]*beta
    se=np.sqrt(np.maximum(0,(residual*residual).sum(0)/df/denominator))
    t=np.divide(beta,se,out=np.zeros_like(beta),where=se>1e-12)
    t[(se<=1e-12)&(np.abs(beta)>1e-12)] = np.inf
    p=2*stats.t.sf(np.abs(t),df)
    return beta,p,df


def top_signature(beta,q,limit):
    candidates=np.flatnonzero(q<=.05)
    return candidates[np.argsort(-np.abs(beta[candidates]),kind='stable')[:limit]]


def pathway_matrix(raw, nodes, gene_ids):
    n=nodes.set_index('index');canonical=nodes.groupby(['type','source','id'])['index'].min()
    pid=nodes[nodes.type.eq('pathway')].sort_values('index').drop_duplicates(['source','id'])
    gene_lookup={str(g):i for i,g in enumerate(gene_ids)}
    path_lookup={str(g):i for i,g in enumerate(pid.id)}
    pairs=set();target_pairs=set()
    use=['x_index','y_index','x_type','y_type','relation','display_relation']
    for chunk in pd.read_csv(raw,usecols=use,chunksize=250000):
        for relation in ['pathway_protein','drug_protein']:
            part=chunk[chunk.relation.eq(relation)]
            if relation=='drug_protein':part=part[part.display_relation.eq('target')]
            for x,y,xt,yt in part[['x_index','y_index','x_type','y_type']].itertuples(index=False,name=None):
                if xt=='gene/protein':g,other=x,y
                elif yt=='gene/protein':g,other=y,x
                else:raise ValueError('Unexpected relation endpoints')
                gene=str(n.loc[g,'id']);other_id=str(n.loc[other,'id'])
                if gene not in gene_lookup:continue
                if relation=='pathway_protein':pairs.add((gene_lookup[gene],path_lookup[other_id]))
                else:target_pairs.add((gene,other_id))
    a,b=np.array(sorted(pairs)).T
    matrix=sparse.csr_matrix((np.ones(len(a)),(a,b)),shape=(len(gene_ids),len(pid)))
    count=np.asarray(matrix.sum(0)).ravel()
    kept=count>=5
    matrix=matrix[:,kept]@sparse.diags(1/count[kept])
    return matrix,pid.iloc[np.flatnonzero(kept)].reset_index(drop=True),target_pairs,count


def model_metrics(name,pred,y,metadata):
    result={'model':name,'accuracy':accuracy_score(y,pred),'balanced_accuracy':balanced_accuracy_score(y,pred),'macro_f1':f1_score(y,pred,labels=CMS,average='macro',zero_division=0)}
    for cls in CMS:
        mask=y==cls;result[cls+'_recall']=float(np.mean(pred[mask]==y[mask]))
    return result


def markdown(frame):
    def fmt(x):return f'{x:.4f}' if isinstance(x,(float,np.floating)) else str(x)
    return '\n'.join(['| '+' | '.join(frame.columns)+' |','| '+' | '.join(['---']*len(frame.columns))+' |',*['| '+' | '.join(fmt(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)]])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'data/CRCSC_data')
    parser.add_argument('--annotation',type=Path,default=ROOT/'data/external/geo/GPL570.annot.gz')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--signature-size',type=int,default=100)
    args=parser.parse_args()
    if args.signature_size<1:parser.error('Signature size must be positive')
    out=args.output or ROOT/'reports/crcsc'/datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    out.mkdir(parents=True,exist_ok=False)
    sources=[*sorted(args.input.iterdir()),args.annotation,ROOT/'data/raw/primekg/kg.csv',ROOT/'data/processed/primekg/nodes.csv']
    hashes={str(p.relative_to(ROOT)):sha(p) for p in sources if p.is_file()}
    print(f'Output: {out}',flush=True)
    annotation,annotation_meta=read_annotation(args.annotation)
    expression,audit,columns=load_expression(args.input/'GSE39582_frma_expression.tsv',annotation)
    audit.to_csv(out/'probe_mapping_audit.csv.gz',index=False,compression='gzip')
    clinical=pd.read_csv(args.input/'FRENCH_clinical.tsv',sep='\t',keep_default_na=False)
    labels=pd.read_csv(args.input/'cms_labels_public_all.txt',sep='\t',keep_default_na=False)
    labels=labels[labels.dataset.eq('gse39582')]
    assert clinical.id.is_unique and labels['sample'].is_unique
    samples=clinical.merge(labels[['sample','CMS_network','CMS_RFclassifier',FINAL]],left_on='id',right_on='sample',validate='1:1')
    assert set(expression.columns)==set(samples.id)
    samples=samples.set_index('id').loc[expression.columns].rename_axis('id').reset_index()
    samples['CMS']=samples[FINAL]
    samples.to_csv(out/'sample_audit.csv',index=False)
    nodes=pd.read_csv(ROOT/'data/processed/primekg/nodes.csv',dtype={'id':str},keep_default_na=False)
    genes=nodes[nodes.type.eq('gene/protein')].sort_values('index').drop_duplicates(['source','id'])
    assert genes.id.is_unique
    map_audit=pd.DataFrame({'ncbi_id':expression.index.astype(str)}).merge(genes[['index','id','name']],left_on='ncbi_id',right_on='id',how='left',validate='1:1')
    map_audit.to_csv(out/'gene_mapping_audit.csv',index=False)
    overlapping=map_audit.loc[map_audit['index'].notna(),'ncbi_id']
    expression=expression.loc[overlapping]
    trainmask=samples.dataset.eq('discovery')&samples.CMS.isin(CMS)
    testmask=samples.dataset.eq('validation')&samples.CMS.isin(CMS)
    train_meta=samples[trainmask].reset_index(drop=True);test_meta=samples[testmask].reset_index(drop=True)
    assert set(train_meta.id).isdisjoint(test_meta.id)
    raw=expression.to_numpy(dtype=np.float64).T
    train,test,variable=standardized(raw[trainmask],raw[testmask])
    gene_ids=expression.index.to_numpy()[variable]
    ytrain=train_meta.CMS.to_numpy();ytest=test_meta.CMS.to_numpy()
    counttable=pd.crosstab(samples.CMS,samples.dataset).reindex(['CMS1','CMS2','CMS3','CMS4','NOLBL'],fill_value=0)
    counttable.to_csv(out/'cohort_counts.csv')
    print(f'Unique mapped/variable genes: {len(gene_ids)}; discovery {len(train)}, validation {len(test)}',flush=True)
    pathway,pmeta,target_pairs,path_counts=pathway_matrix(ROOT/'data/raw/primekg/kg.csv',nodes,gene_ids)
    ptr=train@pathway;pte=test@pathway
    pmeta.to_csv(out/'pathway_inventory.csv',index=False)
    ztrain=covariates(train_meta);ztest=covariates(test_meta)
    gene_tables=[];path_tables=[];replication=[];signatures=[];weights=np.zeros((len(gene_ids),4));targetrows=[]
    for i,subtype in enumerate(CMS):
        discovery_beta,discovery_p,df=adjusted_contrast(train,ytrain==subtype,ztrain)
        validation_beta,validation_p,vdf=adjusted_contrast(test,ytest==subtype,ztest)
        q=bh(discovery_p);valq=bh(validation_p)
        selected=top_signature(discovery_beta,q,args.signature_size)
        weights[selected,i]=discovery_beta[selected]/max(np.abs(discovery_beta[selected]).sum(),1e-12)
        gf=pd.DataFrame({'subtype':subtype,'ncbi_id':gene_ids,'discovery_adjusted_effect':discovery_beta,'discovery_p':discovery_p,'discovery_q':q,'validation_adjusted_effect':validation_beta,'validation_p':validation_p,'validation_q':valq})
        gf=gf.merge(genes[['id','index','name']],left_on='ncbi_id',right_on='id',how='left',validate='1:1')
        gene_tables.append(gf)
        chosen=gf.iloc[selected].copy();chosen['profile_weight']=weights[selected,i];signatures.append(chosen)
        pathway_beta,pp,pdf=adjusted_contrast(ptr,ytrain==subtype,ztrain)
        vp,vpp,vpdf=adjusted_contrast(pte,ytest==subtype,ztest)
        pq=bh(pp);vqq=bh(vpp)
        pf=pmeta[['id','name']].copy();pf['subtype']=subtype;pf['discovery_adjusted_effect']=pathway_beta;pf['discovery_p']=pp;pf['discovery_q']=pq;pf['validation_adjusted_effect']=vp;pf['validation_p']=vpp;pf['validation_q']=vqq;path_tables.append(pf)
        selected_paths=top_signature(pathway_beta,pq,args.signature_size)
        for level,beta,vbeta,qs,vqs,selection in [('gene',discovery_beta,validation_beta,q,valq,selected),('pathway',pathway_beta,vp,pq,vqq,selected_paths)]:
            correlation=stats.spearmanr(beta,vbeta).statistic
            replicate=(np.sign(beta[selection])==np.sign(vbeta[selection]))
            replication.append({'subtype':subtype,'level':level,'tested_features':len(beta),'discovery_q05_features':int((qs<=.05).sum()),'selected_features':len(selection),'all_effect_spearman':float(correlation),'selected_sign_agreement':float(replicate.mean()) if len(selection) else None,'selected_validation_q05_same_direction':int(((vqs[selection]<=.05)&replicate).sum()),'discovery_residual_df':df if level=='gene' else pdf,'validation_residual_df':vdf if level=='gene' else vpdf})
        selected_genes=set(gene_ids[selected]);connected={(g,d) for g,d in target_pairs if g in selected_genes}
        targetrows.append({'subtype':subtype,'selected_genes':len(selected),'selected_genes_with_targets':len({g for g,d in connected}),'unique_connected_drugs':len({d for g,d in connected}),'interpretation':'connectivity only; no action direction, vulnerability or efficacy inferred'})
    pd.concat(gene_tables).to_csv(out/'gene_profiles.csv.gz',index=False,compression='gzip')
    pd.concat(path_tables).to_csv(out/'pathway_profiles.csv.gz',index=False,compression='gzip')
    pd.concat(signatures).to_csv(out/'selected_gene_signatures.csv',index=False)
    pd.DataFrame(replication).to_csv(out/'profile_replication.csv',index=False)
    pd.DataFrame(targetrows).to_csv(out/'signature_drug_connectivity.csv',index=False)
    overlaps=[]
    for i,a in enumerate(CMS):
        for j,b in enumerate(CMS):
            if j<=i:continue
            ga=set(gene_ids[weights[:,i]!=0]);gb=set(gene_ids[weights[:,j]!=0]);union=len(ga|gb)
            overlaps.append({'subtype_a':a,'subtype_b':b,'shared_genes':len(ga&gb),'jaccard':len(ga&gb)/union if union else None})
    pd.DataFrame(overlaps).to_csv(out/'signature_overlap.csv',index=False)
    # Prespecified diagnostics: no test-guided parameter or gene selection.
    results=[];predictions=test_meta[['id','CMS','batch','stage']].copy()
    def assess(name,pred):
        results.append(model_metrics(name,pred,ytest,test_meta));predictions[name]=pred
        pd.DataFrame(confusion_matrix(ytest,pred,labels=CMS),index=CMS,columns=CMS).to_csv(out/f'confusion_{name}.csv')
    assess('majority',np.repeat(pd.Series(ytrain).value_counts().index[0],len(ytest)))
    batch=make_pipeline(OneHotEncoder(handle_unknown='ignore'),LogisticRegression(C=1,max_iter=2000,class_weight='balanced',random_state=42))
    batch.fit(train_meta[['batch']],ytrain);assess('batch_only',batch.predict(test_meta[['batch']]))
    def classify(name,xtrain,xtest):
        model=make_pipeline(StandardScaler(),LogisticRegression(C=1,max_iter=2000,class_weight='balanced',random_state=42))
        model.fit(xtrain,ytrain);assess(name,model.predict(xtest))
    # Variance ranking is estimated from discovery raw gene values only.
    variable_raw=raw[trainmask][:,variable]
    top=np.argsort(-variable_raw.var(0),kind='stable')[:min(2000,len(gene_ids))]
    classify('expression_logistic',train[:,top],test[:,top])
    centroids=np.stack([train[ytrain==c][:,top].mean(0) for c in CMS])
    distance=((test[:,None,top]-centroids[None,:,:])**2).mean(2)
    assess('expression_nearest_centroid',np.array(CMS)[distance.argmin(1)])
    assess('signed_signature_scores',np.array(CMS)[(test@weights).argmax(1)])
    classify('pathway_logistic',ptr,pte)
    for seed in [17,29,43]:
        # Null preserves pathway sizes and membership structure, while permuting
        # which measured gene occupies each position; not per-gene degree fixed.
        perm=np.random.default_rng(seed).permutation(len(gene_ids))
        null=pathway[perm]
        classify('shuffled_pathway_'+str(seed),train@null,test@null)
    result=pd.DataFrame(results);result.to_csv(out/'classification_metrics.csv',index=False);predictions.to_csv(out/'validation_predictions.csv',index=False)
    pd.crosstab(samples.batch,samples.dataset).to_csv(out/'batch_by_split.csv')
    # All validation examples were scored once using prespecified settings.
    # Their selection-independent bootstrap describes this panel, not label truth.
    bootstrap=[]
    rng=np.random.default_rng(42)
    strata=[np.flatnonzero(ytest==c) for c in CMS]
    for repeat in range(1000):
        take=np.concatenate([rng.choice(s,len(s),replace=True) for s in strata])
        for model in predictions.columns[4:]:
            pred=predictions[model].to_numpy()[take]
            baseline=predictions.expression_logistic.to_numpy()[take]
            bootstrap.append({'repeat':repeat,'model':model,'accuracy_delta_vs_expression':float(np.mean(pred==ytest[take])-np.mean(baseline==ytest[take]))})
    bootstrap=pd.DataFrame(bootstrap)
    ci=bootstrap.groupby('model').accuracy_delta_vs_expression.quantile([.025,.5,.975]).unstack();ci.columns=['lower','median','upper'];ci.to_csv(out/'paired_bootstrap_accuracy_delta.csv')
    manifest={'experiment':'CMS molecular profile reproducibility, not drug prediction','annotation_url':ANNOTATION_URL,'annotation_metadata':annotation_meta,'settings':{'signature_size':args.signature_size,'expression_genes':2000,'logistic_C':1,'null_seeds':[17,29,43]},'expression_probe_count':len(audit),'samples':len(samples),'single_id_probes':int(audit.mapping_status.eq('single_ncbi_id').sum()),'unique_annotated_genes':len(map_audit),'primekg_mapped_genes':len(overlapping),'variable_mapped_genes':len(gene_ids),'primekg_mapping_pct':100*len(overlapping)/len(map_audit),'pathways_min5_genes':len(pmeta),'mapping_status':audit.mapping_status.value_counts().to_dict(),'input_sha256':hashes,'input_hashes_unchanged':all(sha(ROOT/p)==h for p,h in hashes.items()),'code_sha256':sha(Path(__file__)),'versions':{'pandas':pd.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__}}
    assert manifest['input_hashes_unchanged']
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    report=f'''# CRCSC CMS molecular-profile experiment

## Scope and outcome

This new experiment uses the uploaded GSE39582 expression/clinical/CMS files. It evaluates whether **discovery-derived gene and pathway profiles reproduce in the supplied validation group** and whether known CMS assignments can be reproduced with simple models. It does not train a GNN, merge graphs, rank drugs, or validate drug–subtype links. The existing CRC agent remains unchanged.

## Cohort and identifiers

All {len(samples)} expression samples join one-to-one to clinical/CMS records by GSM ID. Expression has {len(audit):,} unique Affymetrix probes with finite numeric values. The final combined CMS column is used; NOLBL samples are excluded from supervised profiles/classification, not forced into a class. Network-only and RF-only labels remain in the sample audit.

{markdown(counttable.reset_index())}

Probe annotation: [official GPL570 annotation]({ANNOTATION_URL}); metadata dated Aug 09 2016. Accept only a single distinct numeric NCBI ID, excluding ambiguous multi-gene probes and unmapped/control probes. Aggregate multiple accepted probes per gene by within-sample median. The snapshot is old and mappings are provisional; no alias/name guesses or silent remapping of retired IDs are used. Of {len(map_audit):,} annotated unique genes, **{len(overlapping):,} map by exact NCBI ID to PrimeKG ({manifest['primekg_mapping_pct']:.2f}%)**; {len(gene_ids):,} have nonzero discovery variance. Exact-ID PrimeKG node aliases are canonicalized. The expression header omits the probe column label; the reader explicitly handles the extra first field.

## Discovery-only profile construction

Discovery gene means/standard deviations are used to standardize both groups. For each CMS versus the other three, OLS estimates the expression contrast controlling batch, age, gender, stage, and tumor location, separately within each cohort. BH correction is per subtype and feature level. Up to {args.signature_size} discovery q<=.05 genes, ranked by absolute adjusted effect, form a signed, L1-normalized profile. No validation outcome selects features, weights, or thresholds. Standardization reflects relative cohort expression, **not tumor-versus-normal expression**.

PrimeKG `pathway_protein` pairs define Reactome membership. Pathways with >=5 measured variable genes are retained ({len(pmeta):,}). Features are mean discovery-standardized expression of measured members, not validated activation scores. Overlapping and broad pathways are correlated. Gene permutations preserve pathway sizes and membership structure but do not fix each named gene's degree.

## Profile replication

Spearman compares adjusted effects across all measured features. Sign agreement and validation BH support refer to discovery-selected signatures only; validation labels are used for evaluation, not training. Shared expression-based labeling and processing can make this easier than genuinely independent validation.

{markdown(pd.DataFrame(replication))}

## Prespecified subtype-reproduction diagnostics

Fit on labeled discovery samples, evaluate on labeled validation samples. No test-tuned hyperparameters. Gene logistic/nearest-centroid models use the top 2,000 discovery-variance genes. Logistic models have C=1, class-balanced loss, and max_iter=2,000. StandardScaler is fitted on discovery only. Signature classification uses the four signed expression-profile scores; pathway/null models use the same logistic recipe.

{markdown(result)}

The three shuffled-membership controls test whether curated pathway assignments improve reproduction beyond arbitrary expression pooling at the same feature count. Nulls may retain gene signal, and overlapping pathways make this a restricted control. A high subtype score is **not novel subtype discovery**: CMS labels were themselves derived from expression, and this cohort contributed to the published CMS work. The original discovery/validation grouping is useful, but does not make the CMS labels independently assigned from non-expression evidence. A new external cohort is required for stronger claims. Accuracy deltas against expression logistic are described in `paired_bootstrap_accuracy_delta.csv`; within-class patient resampling on one panel does not establish external generalization or undo label circularity.

## Direct drug-target connectivity of signatures

{markdown(pd.DataFrame(targetrows))}

Only PrimeKG `drug_protein` display role `target` is counted; reciprocal rows and aliases are deduplicated by stable IDs. These are connectivity counts, **not candidate-drug rankings**. Increased/decreased expression does not establish tumor dependency, required pharmacological direction, or therapeutic benefit. A final drug–subtype predictor still needs suitable supervision and independent response evaluation.

## What this supports and does not support

The appropriate conclusion is whether these CMS expression profiles are a reproducible, identifier-audited molecular input. Compare the measured profile consistency and classifier/null results before deciding whether pathway structure adds value. Better CMS reproduction alone cannot establish drug prediction, GNN superiority, genotype causality, patient survival benefit, or clinical efficacy. It also cannot establish a new contribution beyond existing CMS classification.

This cohort has selected mutation-status fields, not the full mutation matrix used in ECoGraph. CMS1–4 must remain separate from GI.CIN/GI.HM-indel/GI.HM-SNV/GI.GS; there is no assumed one-to-one crosswalk. Bulk expression contains immune/stromal composition, particularly relevant to mesenchymal/immune signatures. Batch, stage and location can confound effects despite adjustment; OLS is a descriptive model with limited validation group sizes. FDR and bootstrap intervals do not correct provenance or selection bias. The filename suggests fRMA preprocessing, but raw CEL-level preprocessing was not independently reproduced. Four clinical records have stage0; their interpretation warrants upstream review and no assumption that these are matched normal controls is made.

## Next decision

If profiles replicate, retain them as a **new CMS representation benchmark**, not a drop-in proof of better drug inference. A future GNN can integrate these profiles through mapped genes, but should first demonstrate a benefit over expression/pathway baselines on a clearly defined held-out outcome. For subtype–drug links, obtain independent cell-model/experimental response evidence with defensible CMS assignments. Graph indication learning and literature review alone do not verify therapeutic benefit; literature retrieval should audit support, contradiction, and prior knowledge with citations.

## Files and reproduction

Run `python3 analysis/crcsc_experiment.py`; `--output` must not already exist. Dependencies: pandas, NumPy, SciPy, scikit-learn, threadpoolctl. Official annotation is a local prerequisite; no runtime download occurs. Source files are unchanged (SHA-256 audit: **{manifest['input_hashes_unchanged']}**). See `manifest.json` for versions and hashes.

Outputs: `sample_audit.csv`, `cohort_counts.csv`, `probe_mapping_audit.csv.gz`, `gene_mapping_audit.csv`, `gene_profiles.csv.gz`, `selected_gene_signatures.csv`, `pathway_profiles.csv.gz`, `pathway_inventory.csv`, `profile_replication.csv`, `signature_overlap.csv`, `signature_drug_connectivity.csv`, `classification_metrics.csv`, `validation_predictions.csv`, `confusion_*.csv`, `batch_by_split.csv`, `paired_bootstrap_accuracy_delta.csv`, `manifest.json`, `findings.md`.

Primary context: [CMS study](https://www.nature.com/articles/nm.3967), [GSE39582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582). Local file integrity is checked; exact download provenance of the three supplied files and upstream preprocessing remains to be documented by the team.
'''
    (out/'findings.md').write_text(report)
    print(result.to_string(index=False),flush=True)
    print(pd.DataFrame(replication).to_string(index=False),flush=True)
    print('Completed:',out,flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
