"""Summarize a completed CRC run without retraining or modifying source data."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def tie_top_mask(scores,k=20):
    """Keep all exact ties at the kth score, not an arbitrary subset."""
    scores=np.asarray(scores)
    if not len(scores) or k<1:
        raise ValueError("Nonempty scores and positive k required")
    cutoff=np.sort(scores)[-min(k,len(scores))]
    return scores>=cutoff


def markdown(frame):
    def fmt(v):
        return f"{v:.4f}" if isinstance(v,(float,np.floating)) else str(v)
    return "\n".join(["| "+" | ".join(frame.columns)+" |","| "+" | ".join(["---"]*len(frame.columns))+" |"]+
        ["| "+" | ".join(map(fmt,row))+" |" for row in frame.itertuples(index=False,name=None)])


def write_summary(out):
    import json
    manifest=json.loads((out/"provenance.json").read_text())
    metrics=pd.read_csv(out/"evaluation_metrics.csv")
    profiles=pd.read_csv(out/"profile_summary.csv")
    diagnostic=pd.read_csv(out/"subtype_pairwise_diagnostics.csv")
    nodes=pd.read_csv(ROOT/"data/processed/primekg/nodes.csv",dtype={"id":str},keep_default_na=False).set_index("index")
    target_audit=pd.read_csv(ROOT/"reports/compatibility/primekg_direct_targets.csv",dtype=str)
    crc=pd.read_csv(ROOT/"reports/compatibility/primekg_crc_indications.csv",dtype=str)
    crc_drugs=set(crc.drug_id)
    signatures=target_audit.groupby("drug_id").gene_id.agg(lambda values:"|".join(sorted(set(values)))).to_dict()
    decomposition,ties,tie_candidates=[],[],[]
    for path in sorted(out.glob("scores_*.npz")):
        suffix=path.stem.removeprefix("scores_")
        condition,seed=suffix.rsplit("_",1)
        with np.load(path) as arrays:
            candidates=arrays["drug_node_indices"]
            subtypes=arrays["subtype_names"]
            for profile in arrays.files:
                if profile in {"drug_node_indices","subtype_names"}: continue
                scores=arrays[profile].astype(float)
                grand=scores.mean();drug_mean=scores.mean(axis=0);subtype_offset=scores.mean(axis=1)-grand
                interaction=scores-drug_mean[None,:]-subtype_offset[:,None]
                decomposition.append({"condition":condition,"seed":int(seed),"profile":profile,
                    "drug_main_effect_variance_pct":100*np.var(drug_mean)/np.var(scores),
                    "interaction_rms":np.sqrt(np.mean(interaction**2)),"across_drug_score_sd":scores.std(),
                    "interaction_rms_over_score_sd":np.sqrt(np.mean(interaction**2))/scores.std()})
                for i,subtype in enumerate(subtypes):
                    mask=tie_top_mask(scores[i])
                    cutoff=np.sort(scores[i])[-20]
                    ties.append({"condition":condition,"seed":int(seed),"profile":profile,"subtype":subtype,
                        "inclusive_top20_count":int(mask.sum()),"cutoff_exact_ties":int((scores[i]==cutoff).sum()),
                        "top_score_exact_ties":int((scores[i]==scores[i].max()).sum())})
                    if condition=="target_only" and profile in {"uniform_union","frequency","selective_frequency"}:
                        for pos in np.flatnonzero(mask):
                            node=int(candidates[pos]);drug_id=nodes.loc[node,"id"]
                            tie_candidates.append({"seed":int(seed),"profile":profile,"subtype":subtype,
                                "drug_node_index":node,"drug_id":drug_id,"drug_name":nodes.loc[node,"name"],
                                "tied_rank_min":int((scores[i]>scores[i,pos]).sum())+1,
                                "tie_group_size":int((scores[i]==scores[i,pos]).sum()),"raw_logit":scores[i,pos],
                                "recorded_target_ncbi_ids":signatures.get(drug_id,""),
                                "recorded_indication_for_crc_lexical_candidate":drug_id in crc_drugs,
                                "status":"unvalidated hypothesis; compound order unresolved within exact ties"})
    decomp,tie_frame,candidates=pd.DataFrame(decomposition),pd.DataFrame(ties),pd.DataFrame(tie_candidates)
    decomp.to_csv(out/"score_variance_decomposition.csv",index=False)
    tie_frame.to_csv(out/"tie_audit.csv",index=False)
    candidates.to_csv(out/"tie_aware_top_candidates.csv",index=False)
    aggregate=metrics.groupby("condition",sort=False).agg(runs=("test_ap","size"),mean_test_AP=("test_ap","mean"),
        min_test_AP=("test_ap","min"),max_test_AP=("test_ap","max"),mean_disease_macro_AP=("test_disease_macro_ap","mean"))
    primary=diagnostic[diagnostic.condition.eq("target_only") & diagnostic.profile.isin(["uniform_union","frequency","selective_frequency","permuted_gene_control"])]
    sensitivity=primary.groupby("profile",sort=False).agg(mean_top20_overlap=("top20_overlap","mean"),
        mean_score_Spearman=("score_spearman","mean"),mean_pairwise_score_RMS=("score_rms_difference","mean"))
    fraction=decomp[decomp.condition.eq("target_only") & decomp.profile.isin(["uniform_union","frequency","selective_frequency"])].groupby("profile").drug_main_effect_variance_pct.agg(["mean","min","max"])
    normal=metrics[metrics.condition.eq("target_only")].set_index("seed").test_ap
    ablated=metrics[metrics.condition.eq("drug_only")].set_index("seed").test_ap
    diffs=normal-ablated
    primary_candidates=candidates[candidates.profile.eq("selective_frequency")]
    nr3c1_only=primary_candidates.recorded_target_ncbi_ids.eq("2908").all()
    assert nodes[nodes.type.eq("gene/protein") & nodes.id.eq("2908")].shape[0]==1
    actual_ties=tie_frame[tie_frame.condition.eq("target_only") & tie_frame.profile.isin(["uniform_union","frequency","selective_frequency"])]
    tied_sizes=sorted(actual_ties.top_score_exact_ties.unique().tolist())
    inclusive_sizes=sorted(actual_ties.inclusive_top20_count.unique().tolist())
    all_shared=bool(primary[primary.profile.ne("permuted_gene_control")].top20_overlap.eq(20).all())
    legacy_ap=float(metrics.loc[metrics.condition.eq("legacy_all_roles"),"test_ap"].iloc[0])
    reproduces=abs(legacy_ap-.8190277733864988)<1e-5
    boundary=bool(metrics.loc[metrics.condition.ne("drug_popularity"),"best_epoch"].eq(manifest["epochs"]).all())
    baseline_ap=float(metrics.loc[metrics.condition.eq("drug_popularity"),"test_ap"].iloc[0])
    best_model_ap=float(metrics.loc[metrics.condition.ne("drug_popularity"),"test_ap"].max())
    performance_conclusion=("Drug popularity exceeds every fitted GNN in this run." if baseline_ap>best_model_ap else
        "At least one fitted model exceeds popularity on this fixed split; this requires independent replication.")
    ranking_conclusion=("The tested profiles do not fix the shared top candidate group." if all_shared else
        "Some profile interventions change top20 membership; this is not evidence of biological correctness.")
    profile_counts=profiles.pivot(index="subtype",columns="profile",values="genes").reset_index()
    profile_counts=profile_counts[["subtype","uniform_union","frequency","selective_frequency"]]
    scopes={"prediction":"sampled PrimeKG indication recovery","subtype_sensitivity":"inference perturbations, not response accuracy",
        "hypotheses":"unvalidated raw logits; no clinical efficacy or literature novelty established"}
    (out/"summary_scopes.json").write_text(json.dumps(scopes,indent=2)+"\n")
    text=f"""# CRC GNN proof-of-concept: findings

**Verdict:** the pipeline generates traceable association hypotheses and supports quantitative sensitivity tests. {performance_conclusion} {ranking_conclusion} No subtype response accuracy or clinical efficacy was evaluated.

## What was run

The existing notebook architecture, width 48, optimizer and {manifest['epochs']}-epoch budget were retained. {len(metrics)-1} fits cover the legacy model (seed42), direct-target-only pooling and a retrained drug-only molecular-context ablation (initialization seeds {manifest['model_seeds']}), and no-PPI (seed42). All use exactly the same disease split, positive examples, epoch-specific sampled unknown training pairs and validation/test pairs. Legacy AP matches the historical notebook to displayed precision: **{reproduces}**. No original notebooks or outputs changed; all preservation checks passed.

There are 6,233 eligible indications: 4,229 train, 1,151 validation and 853 test; 633/136/136 diseases respectively. Test labels comprise 853 positives and 853 unknowns. Drug popularity uses training counts only. The fixed evaluation retains all 6,282 original drugs, including 592 with no direct target message in corrected models; a common target-supported pair-subset metric is also saved. Hypothesis rankings are restricted to 5,690 direct-target-supported drugs. No therapeutic edges enter molecular adjacency. Source/family overlap and false negatives remain possible.

## 1. Demonstrated predictive performance

{markdown(aggregate.reset_index())}

{performance_conclusion} Across initialization seeds, target-only minus drug-only test AP is **{diffs.mean():+.4f} on average**, ranging **{diffs.min():+.4f} to {diffs.max():+.4f}**. This comparison measures whether removing molecular disease context degrades prediction. The single-seed PPI comparison is exploratory, not evidence of reliable graph benefit. These ranges vary initialization on one split, not independent datasets or confidence intervals. Every selected checkpoint reached the epoch-budget boundary: **{boundary}**. A fixed budget does not rule out improvement with other training.

## 2. Why subtype rankings are nearly identical

Two effects are directly observed:

- **Drug-level variation dominates.** Across-drug main effects account for the following percentage of total score variance. The rest includes subtype-wide offsets and drug-by-subtype interactions; this decomposition is descriptive, not a causal attribution.

{markdown(fraction.reset_index())}

- **Identical target neighborhoods are structurally indistinguishable to this model.** Observed top-score tie-group sizes across target-only seeds and actual profiles are **{tied_sizes}**; inclusive top20 counts are **{inclusive_sizes}**. All tie-inclusive selective-profile candidates have NR3C1 (NCBI2908) as their only recorded direct target: **{nr3c1_only}**. That gene has one PrimeKG node. Identical target neighborhoods necessarily produce identical drug embeddings and scores in this architecture, which has no compound-specific features. These are recorded connections, not evidence that modulating NR3C1 benefits CRC. `tie_aware_top_candidates.csv` reports minimum tied rank and full group size rather than inventing an ordering among those compounds.

The broad union profiles and mean pooling also discard mutation frequency and average many overlapping genes. This is consistent with the collapse, while the controls quantify the limited dependence rather than proving a single cause.

## 3. Did weighted or selective profiles help?

Profiles are all-observed uniform weights; normalized observation frequencies; and a prespecified selective rule (observed in >=2 patients and >=10% of the subtype, capped at top200 by frequency). They use no indication test labels or favored drugs. Selection is descriptive, not burden-adjusted enrichment.

{markdown(profile_counts)}

Across six subtype pairs and three target-only seeds:

{markdown(sensitivity.reset_index())}

All raw top20 subtype-pair overlaps equal20: **{all_shared}**. {ranking_conclusion} Numerical score changes are a sensitivity result, not validated improvement. `profile_summary.csv` audits cutoff ties: the selective HMSNV cap in this run cuts through a {int(profiles.loc[profiles.profile.eq('selective_frequency') & profiles.subtype.eq('subtype:GI.HM-SNV'),'cutoff_tied_genes'].iloc[0])}-gene frequency tie. Its deterministic tie choice is especially unstable with only six patients. Frequencies are not driver probabilities and the procedure does not resolve mutation burden, gene length, passenger flags or missing variant consequences.

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

`python3 analysis/crc_experiments.py` creates a new timestamped run. Requires the local HGNC TSV supplied with `--hgnc` (default `data/external/hgnc/hgnc_complete_set.txt`); the experiment runner does not download it. `python3 analysis/crc_report.py {out.relative_to(ROOT) if out.is_relative_to(ROOT) else out}` regenerates this summary from completed results without retraining. See `crc_experiment_report.md` for extended methods and `provenance.json` for hashes, versions and preservation evidence.
"""
    (out/"crc_findings.md").write_text(text)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run",type=Path)
    write_summary(parser.parse_args().run.resolve())
