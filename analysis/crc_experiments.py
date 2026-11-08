"""Incremental CRC GNN diagnostics, preserving original notebooks and outputs.

Run: python3 analysis/crc_experiments.py
No live inference, clinical efficacy labels, or full graph merge is performed.
"""

import argparse
import hashlib
import json
from copy import deepcopy
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn
import torch
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score
from torch import nn
from torch.nn import functional as F

ROOT = Path(__file__).resolve().parents[1]
LABEL_RELATIONS = ["indication", "off-label use", "contraindication"]
PROFILE_NAMES = ["uniform_union", "frequency", "selective_frequency"]


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def weighted_adjacency(left, right, weights, n_left, n_right):
    """Row-normalized sparse weighted pooling; empty rows stay zero."""
    left, right, weights = np.asarray(left), np.asarray(right), np.asarray(weights, dtype=float)
    if not (len(left) == len(right) == len(weights)):
        raise ValueError("Mismatched sparse edge arrays")
    if np.any(weights < 0) or not np.isfinite(weights).all():
        raise ValueError("Pooling weights must be finite and nonnegative")
    totals = np.bincount(left.astype(int), weights=weights, minlength=n_left)
    normalized = np.divide(weights, totals[left.astype(int)], out=np.zeros_like(weights), where=totals[left.astype(int)] > 0)
    return torch.sparse_coo_tensor(
        torch.tensor(np.stack([left, right]), dtype=torch.long),
        torch.tensor(normalized, dtype=torch.float32), (n_left, n_right),
    ).coalesce()


def make_profiles(matrix, labels, genes, mapped):
    """Observation profiles only: no enrichment/causality claims or outcome tuning."""
    profiles, rows, evidence = {}, [], []
    subtypes = sorted(labels.unique())
    gene_pos = pd.Series(np.arange(len(genes)), index=genes)
    for name in PROFILE_NAMES:
        array = np.zeros((len(subtypes), len(genes)), dtype=np.float32)
        for i, subtype in enumerate(subtypes):
            patients = labels.index[labels.eq(subtype)]
            counts = matrix.loc[patients, mapped.index].sum().astype(int)
            frequencies = counts / len(patients)
            chosen = counts[counts.gt(0)].index
            cutoff_ties = 0
            if name == "selective_frequency":
                # Prespecified descriptive rule, not chosen using prediction outputs.
                threshold = max(2, int(np.ceil(.1 * len(patients))))
                eligible = pd.DataFrame({"count": counts, "frequency": frequencies,
                    "node": mapped.node_index}).query("count >= @threshold")
                eligible = eligible.sort_values(["frequency", "node"], ascending=[False, True], kind="stable")
                if len(eligible) > 200:
                    cutoff_ties = int(eligible.frequency.eq(eligible.frequency.iloc[199]).sum())
                chosen = eligible.head(200).index
            positions = gene_pos.loc[mapped.loc[chosen, "node_index"]].to_numpy()
            weights = np.ones(len(chosen)) if name == "uniform_union" else frequencies.loc[chosen].to_numpy()
            if len(weights) == 0 or weights.sum() == 0:
                raise ValueError(f"Empty profile: {name}, {subtype}")
            array[i, positions] = weights / weights.sum()
            rows.append({"profile": name, "subtype": subtype, "samples": len(patients),
                "genes": len(chosen), "effective_genes": 1 / np.square(array[i]).sum(),
                "cutoff_tied_genes": cutoff_ties})
            for symbol, weight in zip(chosen, weights):
                evidence.append({"profile": name, "subtype": subtype, "symbol": symbol,
                    "gene_node_index": int(mapped.loc[symbol, "node_index"]),
                    "n_observed": int(counts[symbol]), "n_profiled": len(patients),
                    "frequency": frequencies[symbol], "unnormalized_weight": weight,
                    "pooling_weight": weight / weights.sum()})
        profiles[name] = array
    return subtypes, profiles, pd.DataFrame(rows), pd.DataFrame(evidence)


class IndicationGNN(nn.Module):
    """Notebook architecture; switches isolate graph/context information."""
    def __init__(self, n_genes, adj_ppi, adj_drug, adj_disease, use_context=True, use_ppi=True, width=48):
        super().__init__()
        self.adj_ppi, self.adj_drug, self.adj_disease = adj_ppi, adj_drug, adj_disease
        self.use_context, self.use_ppi = use_context, use_ppi
        self.gene = nn.Embedding(n_genes, width)
        self.gene_self = nn.Linear(width, width)
        self.gene_ppi = nn.Linear(width, width, bias=False)
        self.drug_projection = nn.Linear(width, width)
        self.disease_projection = nn.Linear(width, width)
        self.decoder = nn.Sequential(nn.Linear(3 * width, width), nn.ReLU(), nn.Linear(width, 1))

    def encode(self):
        initial = self.gene.weight
        neighbors = torch.sparse.mm(self.adj_ppi, initial) if self.use_ppi else torch.zeros_like(initial)
        gene_vectors = F.relu(self.gene_self(initial) + self.gene_ppi(neighbors))
        drug_vectors = self.drug_projection(torch.sparse.mm(self.adj_drug, gene_vectors))
        pooled = torch.sparse.mm(self.adj_disease, gene_vectors)
        if not self.use_context:
            pooled = torch.zeros_like(pooled)
        return self.disease_projection(pooled), drug_vectors, gene_vectors

    def score(self, encoded, pairs):
        ids = torch.tensor(pairs, dtype=torch.long)
        return self.score_vectors(encoded[0][ids[:, 0]], encoded[1][ids[:, 1]])

    def score_vectors(self, disease, drug):
        return self.decoder(torch.cat([disease, drug, disease * drug], dim=1)).squeeze(1)


def load_inputs(hgnc_path, out):
    # Restore the notebook's endpoint-encounter order to preserve gene embedding ordering.
    metadata = pd.read_csv(ROOT / "data/processed/primekg/nodes.csv", dtype={"id": str}, keep_default_na=False).set_index("index")
    chunks, orders, raw_rows = [], [], 0
    selected = ["protein_protein", "drug_protein", "disease_protein", *LABEL_RELATIONS]
    for chunk in pd.read_csv(ROOT / "data/raw/primekg/kg.csv",
            usecols=["x_index", "y_index", "relation", "display_relation"], chunksize=250_000, keep_default_na=False):
        raw_rows += len(chunk)
        orders.append(pd.concat([chunk.x_index, chunk.y_index]).drop_duplicates())
        chunks.append(chunk[chunk.relation.isin(selected)].copy())
    order = pd.concat(orders).drop_duplicates()
    assert set(order) == set(metadata.index)
    metadata = metadata.loc[order]
    edges = pd.concat(chunks, ignore_index=True)
    kind = metadata.type

    def pairs(relations, left_type, right_type, display=None):
        rows = edges[edges.relation.isin(relations)]
        if display:
            rows = rows[rows.display_relation.eq(display)]
        forward = rows.x_index.map(kind).eq(left_type)
        result = pd.DataFrame({"left": rows.x_index.where(forward, rows.y_index),
                               "right": rows.y_index.where(forward, rows.x_index)}).drop_duplicates()
        assert result.left.map(kind).eq(left_type).all() and result.right.map(kind).eq(right_type).all()
        return result

    ppi = pairs(["protein_protein"], "gene/protein", "gene/protein")
    all_roles = pairs(["drug_protein"], "drug", "gene/protein")
    direct = pairs(["drug_protein"], "drug", "gene/protein", "target")
    gd = pairs(["disease_protein"], "disease", "gene/protein")
    known = pairs(["indication"], "disease", "drug")
    excluded = pairs(LABEL_RELATIONS, "disease", "drug")
    genes = pd.Index(metadata.index[metadata.type.eq("gene/protein")])
    drugs = pd.Index(sorted(all_roles.left.unique()))

    eco = ROOT / "ECoGraph/data"
    mutation = pd.read_csv(eco / "coad_mut_data.csv", keep_default_na=False)
    barcode, patient_id = mutation.pop("Tumor_Sample_Barcode"), mutation.pop("patient_id")
    assert patient_id.is_unique and barcode.str[:12].eq(patient_id).all()
    assert barcode.str.match(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01[A-Z0-9]-").all()
    assert mutation.columns.is_unique and mutation.isin([0, 1]).all().all()
    mutation = mutation.astype("uint8").set_axis(patient_id, axis=0)
    assignment = pd.read_csv(eco / "neptune_csv/edges_has_subtype.csv", keep_default_na=False)
    assert assignment["~label"].eq("HAS_SUBTYPE").all()
    labels = pd.Series(assignment["~to"].values, index=assignment["~from"].str.removeprefix("sample:"))
    assert labels.index.is_unique and labels.index.isin(mutation.index).all()
    roster = pd.read_csv(eco / "neptune_csv/vertices_subtype.csv")["~id"]
    assert set(labels) == set(roster) and len(roster) == 4
    hgnc = pd.read_csv(hgnc_path, sep="\t", usecols=["symbol", "entrez_id", "status", "hgnc_id"], dtype=str, keep_default_na=False)
    approved = hgnc.loc[hgnc.status.eq("Approved"), ["symbol", "entrez_id", "hgnc_id"]].drop_duplicates()
    approved = approved[~approved.symbol.duplicated(keep=False)]
    kg_gene = metadata[metadata.type.eq("gene/protein") & metadata.source.eq("NCBI")].reset_index()[["index", "id"]]
    kg_gene = kg_gene[~kg_gene.id.duplicated(keep=False)].rename(columns={"index": "node_index"})
    mapping = pd.DataFrame({"symbol": mutation.columns}).merge(approved, on="symbol", how="left", validate="1:1")
    mapping = mapping.merge(kg_gene, left_on="entrez_id", right_on="id", how="left", validate="m:1")
    mapping["mapped"] = mapping.node_index.notna()
    duplicate = mapping.loc[mapping.mapped, "node_index"].duplicated(keep=False)
    mapping.loc[duplicate.index[duplicate], "mapped"] = False
    mapping.to_csv(out / "gene_mapping_audit.csv", index=False)
    mapped = mapping[mapping.mapped].copy().set_index("symbol")
    mapped["node_index"] = mapped.node_index.astype(int)
    subtypes, profiles, profile_summary, evidence = make_profiles(mutation, labels, genes, mapped)
    subtype_indices = np.arange(int(metadata.index.max()) + 1, int(metadata.index.max()) + 5)
    diseases = pd.Index(sorted(set(gd.left) | set(subtype_indices)))
    gpos = pd.Series(np.arange(len(genes)), index=genes)
    dpos = pd.Series(np.arange(len(diseases)), index=diseases)
    drugpos = pd.Series(np.arange(len(drugs)), index=drugs)
    both = pd.concat([ppi, ppi.rename(columns={"left": "right", "right": "left"})]).drop_duplicates()
    adj_ppi = weighted_adjacency(gpos.loc[both.left], gpos.loc[both.right], np.ones(len(both)), len(genes), len(genes))
    adj_drugs = {}
    for name, rows in [("all_roles", all_roles), ("target", direct)]:
        adj_drugs[name] = weighted_adjacency(drugpos.loc[rows.left], gpos.loc[rows.right], np.ones(len(rows)), len(drugs), len(genes))
    # Keep the original unweighted union in the training graph. Subtypes have no labels,
    # and messages flow gene -> disease only, so their inference profiles cannot affect training.
    sr, sc = np.nonzero(profiles["uniform_union"])
    left = np.r_[dpos.loc[gd.left].values, dpos.loc[subtype_indices].values[sr]]
    right = np.r_[gpos.loc[gd.right].values, sc]
    adj_disease = weighted_adjacency(left, right, np.ones(len(left)), len(diseases), len(genes))
    positives = known[known.left.isin(diseases) & known.right.isin(drugs)].copy()
    target_mask = drugs.isin(direct.left)
    summary = {"raw_rows": raw_rows, "genes": len(genes), "drugs_all_roles": len(drugs),
        "drugs_direct_target": int(target_mask.sum()), "zero_target_drugs_retained_for_matched_evaluation": int((~target_mask).sum()),
        "mapped_matrix_genes": len(mapped), "known_indication_pairs": len(known), "eligible_positive_pairs": len(positives),
        "direct_target_pairs": len(direct), "all_role_drug_gene_pairs": len(all_roles), "typed_patients": len(labels)}
    profile_summary.to_csv(out / "profile_summary.csv", index=False)
    evidence.to_csv(out / "profile_gene_evidence.csv.gz", index=False)
    return locals()


def make_split(data, out, epochs, split_seed=42):
    rng = np.random.default_rng(split_seed)
    eligible = rng.permutation(data["positives"].left.unique())
    size = max(1, round(.15 * len(eligible)))
    identifiers = {"train": eligible[2*size:], "validation": eligible[:size], "test": eligible[size:2*size]}
    diseases, drugs = data["diseases"], data["drugs"]
    positives, negatives, split_rows = {}, {}, []
    for name, ids in identifiers.items():
        rows = data["positives"][data["positives"].left.isin(ids)]
        positives[name] = np.column_stack([diseases.get_indexer(rows.left), drugs.get_indexer(rows.right)])
        split_rows.extend({"split": name, "disease_node_index": int(i)} for i in ids)
    assert set(identifiers["train"]).isdisjoint(identifiers["test"])
    assert set(identifiers["train"]).isdisjoint(identifiers["validation"])
    assert set(identifiers["test"]).isdisjoint(identifiers["validation"])
    forbidden = {d: set(drugs.get_indexer(rows.right)) - {-1} for d, rows in data["excluded"].groupby("left")}
    pools = {diseases.get_loc(d): np.array(sorted(set(range(len(drugs))) - forbidden.get(d, set()))) for d in eligible}
    def sample(rows):
        sampled = []
        for disease, n in zip(*np.unique(rows[:, 0], return_counts=True)):
            if len(pools[disease]) < n:
                raise ValueError("Insufficient unknown pairs for matched sampling")
            sampled.extend((disease, drug) for drug in rng.choice(pools[disease], n, replace=False))
        return np.asarray(sampled, dtype=np.int64)
    negatives["validation"] = sample(positives["validation"])
    negatives["test"] = sample(positives["test"])
    train_negatives = [sample(positives["train"]) for _ in range(epochs)]
    for name in ["validation", "test"]:
        assert not set(map(tuple, positives[name])) & set(map(tuple, negatives[name]))
    pd.DataFrame(split_rows).to_csv(out / "disease_split.csv", index=False)
    np.savez_compressed(out / "evaluation_pairs.npz", **{f"{k}_positive": v for k,v in positives.items()},
                        **{f"{k}_unlabeled": v for k,v in negatives.items()})
    return positives, negatives, train_negatives, identifiers


def pair_metrics(model, positive, negative, target_mask):
    pairs = np.concatenate([positive, negative])
    labels = np.r_[np.ones(len(positive)), np.zeros(len(negative))]
    with torch.no_grad():
        scores = model.score(model.encode(), pairs).numpy()
    supported = target_mask[pairs[:, 1]]
    macro = np.mean([average_precision_score(labels[pairs[:,0] == d], scores[pairs[:,0] == d]) for d in np.unique(pairs[:,0])])
    return {"ap": average_precision_score(labels, scores), "disease_macro_ap": macro,
            "target_supported_ap": average_precision_score(labels[supported], scores[supported])}, scores


def fit(data, split, condition, seed, epochs, out):
    positives, negatives, train_negatives, _ = split
    torch.manual_seed(seed)
    roles = "all_roles" if condition == "legacy_all_roles" else "target"
    model = IndicationGNN(len(data["genes"]), data["adj_ppi"], data["adj_drugs"][roles], data["adj_disease"],
        use_context=condition != "drug_only", use_ppi=condition != "no_ppi")
    optimizer = torch.optim.AdamW(model.parameters(), lr=.003, weight_decay=.001)
    best_ap, best_epoch, best_state = -1., 0, None
    curve = []
    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        pairs = np.concatenate([positives["train"], train_negatives[epoch]])
        labels = np.r_[np.ones(len(positives["train"])), np.zeros(len(train_negatives[epoch]))]
        loss = F.binary_cross_entropy_with_logits(model.score(model.encode(), pairs), torch.tensor(labels, dtype=torch.float32))
        loss.backward()
        optimizer.step()
        model.eval()
        with torch.no_grad():
            pv = np.concatenate([positives["validation"], negatives["validation"]])
            yv = np.r_[np.ones(len(positives["validation"])), np.zeros(len(negatives["validation"]))]
            val_ap = average_precision_score(yv, model.score(model.encode(), pv).numpy())
        curve.append({"condition": condition, "seed": seed, "epoch": epoch+1, "loss": loss.item(), "validation_ap": val_ap})
        if val_ap > best_ap:
            best_ap, best_epoch, best_state = val_ap, epoch+1, deepcopy(model.state_dict())
    model.load_state_dict(best_state)
    model.eval()
    metrics, scores = pair_metrics(model, positives["test"], negatives["test"], data["target_mask"])
    row = {"condition": condition, "seed": seed, "best_epoch": best_epoch, "validation_ap": best_ap,
           "test_ap": metrics["ap"], "test_disease_macro_ap": metrics["disease_macro_ap"],
           "test_target_supported_ap": metrics["target_supported_ap"]}
    print(json.dumps(row), flush=True)
    return model, row, curve, scores


def stable_order(scores, ids):
    """Ties use fixed integer node IDs; differences within ties aren't evidence."""
    return np.lexsort((np.asarray(ids), -np.asarray(scores)))


def rank_diagnostics(scores, drug_ids, subtypes, condition, seed, profile):
    rows = []
    for i,j in combinations(range(len(subtypes)),2):
        topi, topj = stable_order(scores[i], drug_ids)[:20], stable_order(scores[j], drug_ids)[:20]
        a,b = set(topi),set(topj)
        rows.append({"condition": condition, "seed": seed, "profile": profile,
            "subtype_a": subtypes[i], "subtype_b": subtypes[j],
            "top20_overlap": len(a&b), "top20_jaccard": len(a&b)/len(a|b),
            "score_spearman": float(spearmanr(scores[i], scores[j]).statistic),
            "score_rms_difference": float(np.sqrt(np.mean((scores[i]-scores[j])**2)))})
    return rows


def subtype_diagnostics(model, data, condition, seed, out):
    with torch.no_grad():
        encoded = model.encode()
        gene, drug = encoded[2], encoded[1]
        profiles = dict(data["profiles"])
        primary = profiles["selective_frequency"]
        profiles["identical_pooled_control"] = np.repeat(primary.mean(axis=0, keepdims=True), 4, axis=0)
        profiles["zero_context_control"] = np.zeros_like(primary)
        permutation = np.random.default_rng(2026).permutation(primary.shape[1])
        profiles["permuted_gene_control"] = primary[:, permutation]
        profiles["swapped_subtype_control"] = primary[[1,0,3,2]]
        candidates = np.flatnonzero(data["target_mask"])
        candidate_nodes = data["drugs"][candidates]
        score_arrays, vector_arrays, diagnostic_rows, ranks, perturbations = {}, {}, [], [], []
        metadata = data["metadata"]
        for name, weights in profiles.items():
            pooled = torch.tensor(weights) @ gene
            if not model.use_context:
                pooled = torch.zeros_like(pooled)
            vectors = model.disease_projection(pooled)
            scores = np.stack([model.score_vectors(v.expand(len(candidates), -1), drug[candidates]).numpy() for v in vectors])
            score_arrays[name], vector_arrays[name] = scores, vectors.numpy()
            diagnostic_rows.extend(rank_diagnostics(scores, candidate_nodes, data["subtypes"], condition, seed, name))
            # Residual above the same drug's mean score across the four subtypes.
            centered = scores - scores.mean(axis=0,keepdims=True)
            for i, subtype in enumerate(data["subtypes"]):
                for ranking_type, values in [("raw", scores[i]), ("subtype_contrast", centered[i])]:
                    order = stable_order(values, candidate_nodes)[:20]
                    for rank, position in enumerate(order,1):
                        node = candidate_nodes[position]
                        ranks.append({"condition": condition, "seed": seed, "profile": name,
                            "ranking_type": ranking_type, "subtype": subtype, "rank": rank,
                            "drug_node_index": int(node), "drug_id": metadata.loc[node,"id"],
                            "drug_name": metadata.loc[node,"name"], "ranking_score": float(values[position]),
                            "raw_logit": float(scores[i,position]),
                            "interpretation": ("unvalidated association hypothesis; not efficacy or calibrated probability"
                                if name in PROFILE_NAMES and ranking_type == "raw" and condition == "target_only"
                                else "diagnostic control/contrast only; not a candidate recommendation")})
        for name in ["uniform_union", "frequency", "identical_pooled_control", "zero_context_control", "permuted_gene_control"]:
            for i, subtype in enumerate(data["subtypes"]):
                a,b=score_arrays["selective_frequency"][i],score_arrays[name][i]
                overlap=len(set(stable_order(a,candidate_nodes)[:20]) & set(stable_order(b,candidate_nodes)[:20]))
                perturbations.append({"condition":condition,"seed":seed,"subtype":subtype,
                    "reference_profile":"selective_frequency","comparison_profile":name,
                    "top20_overlap":overlap,"score_spearman":float(spearmanr(a,b).statistic),
                    "rms_score_change":float(np.sqrt(np.mean((a-b)**2)))})
        assert np.allclose(score_arrays["swapped_subtype_control"], score_arrays["selective_frequency"][[1,0,3,2]], atol=1e-6)
        assert np.allclose(score_arrays["identical_pooled_control"], score_arrays["identical_pooled_control"][[1,0,3,2]], atol=1e-6)
        embedding_rows=[]
        for name,vectors in vector_arrays.items():
            scores=score_arrays[name]
            for i,j in combinations(range(4),2):
                cosine=float(np.dot(vectors[i],vectors[j])/(np.linalg.norm(vectors[i])*np.linalg.norm(vectors[j])+1e-12))
                embedding_rows.append({"condition":condition,"seed":seed,"profile":name,
                    "subtype_a":data["subtypes"][i],"subtype_b":data["subtypes"][j],
                    "embedding_cosine":cosine,
                    "subtype_score_rms_over_drug_score_sd":float(np.sqrt(np.mean((scores[i]-scores[j])**2))/(scores.std()+1e-12))})
        np.savez_compressed(out/f"scores_{condition}_{seed}.npz",drug_node_indices=candidate_nodes.to_numpy(),
            subtype_names=np.asarray(data["subtypes"]),**score_arrays)
    return diagnostic_rows,ranks,perturbations,embedding_rows


def md_table(frame):
    def format_value(v):
        if isinstance(v,(float,np.floating)): return f"{v:.4f}"
        return str(v)
    return "\n".join(["| "+" | ".join(frame.columns)+" |", "| "+" | ".join(["---"]*len(frame.columns))+" |"]+
        ["| "+" | ".join(map(format_value,row))+" |" for row in frame.itertuples(index=False,name=None)])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hgnc",type=Path,default=ROOT/"data/external/hgnc/hgnc_complete_set.txt")
    parser.add_argument("--seeds",type=int,nargs="+",default=[42,43,44])
    parser.add_argument("--epochs",type=int,default=40)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    if args.epochs<1 or not args.seeds: parser.error("Positive epochs and at least one seed required")
    if not args.hgnc.is_file(): parser.error("Supply the HGNC TSV used by the notebook via --hgnc")
    run_id=datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    out=args.output or ROOT/"reports/crc_experiments"/run_id
    if out.exists(): parser.error("Output directory must be new to preserve previous results")
    protected=[p for base in ["notebooks","reports","artifacts"] for p in (ROOT/base).rglob("*") if p.is_file()]
    protected += [ROOT/"data/raw/primekg/kg.csv",ROOT/"data/processed/primekg/nodes.csv",args.hgnc,
                  ROOT/"ECoGraph/data/coad_mut_data.csv",ROOT/"ECoGraph/data/neptune_csv/edges_has_subtype.csv",
                  ROOT/"ECoGraph/data/neptune_csv/vertices_subtype.csv"]
    original_hashes={str(p):fingerprint(p) for p in protected}
    out.mkdir(parents=True)
    torch.set_num_threads(2)
    torch.use_deterministic_algorithms(True)
    data=load_inputs(args.hgnc,out)
    print("INPUTS",json.dumps(data["summary"]),flush=True)
    split=make_split(data,out,args.epochs)
    positive,negative,_,identifiers=split
    split_summary={name:{"diseases":len(ids),"positive_pairs":len(positive[name])} for name,ids in identifiers.items()}
    pairs=np.concatenate([positive["test"],negative["test"]])
    y=np.r_[np.ones(len(positive["test"])),np.zeros(len(negative["test"]))]
    popularity=np.bincount(positive["train"][:,1],minlength=len(data["drugs"]))
    supported=data["target_mask"][pairs[:,1]]
    baseline={"condition":"drug_popularity","seed":42,"best_epoch":0,"validation_ap":np.nan,
        "test_ap":average_precision_score(y,popularity[pairs[:,1]]),
        "test_disease_macro_ap":np.mean([average_precision_score(y[pairs[:,0]==d],popularity[pairs[pairs[:,0]==d,1]]) for d in np.unique(pairs[:,0])]),
        "test_target_supported_ap":average_precision_score(y[supported],popularity[pairs[supported,1]])}
    metrics,curves,diagnostics,rankings,perturbations,embeddings=[],[],[],[],[],[]
    eval_predictions={"labels":y,"disease_positions":pairs[:,0],"drug_positions":pairs[:,1],"drug_popularity":popularity[pairs[:,1]]}
    schedule=[("legacy_all_roles",42)]+[(condition,seed) for seed in args.seeds for condition in ["target_only","drug_only"]]+[("no_ppi",42)]
    for condition,seed in schedule:
        model,row,curve,scores=fit(data,split,condition,seed,args.epochs,out)
        metrics.append(row);curves.extend(curve);eval_predictions[f"{condition}_{seed}"]=scores
        a,b,c,d=subtype_diagnostics(model,data,condition,seed,out)
        diagnostics.extend(a);rankings.extend(b);perturbations.extend(c);embeddings.extend(d)
        pd.DataFrame([baseline,*metrics]).to_csv(out/"evaluation_metrics.csv",index=False)
        pd.DataFrame(diagnostics).to_csv(out/"subtype_pairwise_diagnostics.csv",index=False)
        pd.DataFrame(rankings).to_csv(out/"hypothesis_rankings.csv",index=False)
        pd.DataFrame(perturbations).to_csv(out/"profile_perturbations.csv",index=False)
        pd.DataFrame(embeddings).to_csv(out/"embedding_diagnostics.csv",index=False)
        pd.DataFrame(curves).to_csv(out/"training_curves.csv",index=False)
    np.savez_compressed(out/"test_predictions.npz",**eval_predictions)
    preservation={str(p):fingerprint(p)==original_hashes[str(p)] for p in protected}
    assert all(preservation.values()),"Protected input/output changed"
    manifest={"input_counts":data["summary"],"split":split_summary,"split_seed":42,
        "model_seeds":args.seeds,"epochs":args.epochs,"training_architecture":"Notebook width48, one PPI layer, mean pooling, same decoder/AdamW",
        "profile_selection":"All union; normalized frequency; count >= max(2,ceil(.1*n)) then top200 frequency, fixed gene-index ties",
        "preservation_checks":preservation,"protected_sha256":original_hashes,
        "versions":{"torch":torch.__version__,"numpy":np.__version__,"pandas":pd.__version__,"sklearn":sklearn.__version__,"scipy":scipy.__version__},
        "clinical_validation":False,"subtype_labels_used_for_indication_training":False,
        "candidate_policy":"All legacy drugs for fixed held-out evaluation; direct-target-supported drugs only for subtype rankings"}
    (out/"provenance.json").write_text(json.dumps(manifest,indent=2)+"\n")
    build_report(out,data,manifest)
    print("RESULTS",out,flush=True)


def build_report(out,data,manifest):
    metrics=pd.read_csv(out/"evaluation_metrics.csv")
    diagnostic=pd.read_csv(out/"subtype_pairwise_diagnostics.csv")
    perturb=pd.read_csv(out/"profile_perturbations.csv")
    aggregate=metrics.groupby("condition",sort=False).agg(runs=("test_ap","size"),mean_test_ap=("test_ap","mean"),
        min_test_ap=("test_ap","min"),max_test_ap=("test_ap","max"),disease_macro_ap=("test_disease_macro_ap","mean"))
    sensitivity=diagnostic[diagnostic.condition.eq("target_only")].groupby("profile",sort=False).agg(
        mean_top20_overlap=("top20_overlap","mean"),mean_score_spearman=("score_spearman","mean"),mean_rms_difference=("score_rms_difference","mean"))
    changes=perturb[perturb.condition.eq("target_only")].groupby("comparison_profile",sort=False).agg(
        mean_top20_overlap=("top20_overlap","mean"),mean_score_spearman=("score_spearman","mean"),mean_rms_change=("rms_score_change","mean"))
    report=f"""# CRC hypothesis-generation experiments

## Scope and minimal changes

Separate scripted experiments reuse the notebook's architecture (width 48, one PPI propagation step, drug/disease gene pooling and the same decoder), 40-epoch default, AdamW settings and disease-held-out protocol. No original notebooks or saved outputs were modified. No full merged graph artifact was available locally, so compact model inputs were reconstructed from the unchanged raw PrimeKG and ECoGraph matrix using the notebook's HGNC mapping rule. The HGNC reference is a new downloaded snapshot, fingerprinted in provenance; it need not equal the teammate's historical snapshot.

The changes are weighted/selective **inference-time** subtype profiles, direct-target filtering as a separately identified model condition, and context/PPI ablations. Subtypes have no supervised indication examples. In this feed-forward model, changing only their outgoing aggregation rows cannot alter training of labeled diseases, so retraining each profile would add no experimental information. The existing model is retained rather than replaced.

## Inputs and fixed evaluation

Input counts: `{json.dumps(data['summary'])}`.

Disease split (seed 42): `{json.dumps(manifest['split'])}`. Every condition uses the exact same train positives, per-epoch sampled training pairs, validation pairs and test pairs. Model seeds {manifest['model_seeds']} vary initialization, not split or patient sampling. Validation AP alone chooses an epoch. Legacy and no-PPI are single-seed diagnostics. Original target/reverse rows are canonicalized; indication/off-label/contraindication labels are excluded from molecular adjacency and from sampled unknown pools. Unknown pairs are not biological negatives. Only exact node identities are split; source overlap and disease-family equivalence remain unresolved.

The candidate universe for disease evaluation deliberately retains **all legacy drug-protein-connected drugs** so correcting target semantics does not silently change the task. Target-only models therefore have {data['summary']['zero_target_drugs_retained_for_matched_evaluation']} drugs with zero target messages; `test_target_supported_ap` separately reports the same pair subset supported by actual target edges. CRC hypothesis rankings use only {data['summary']['drugs_direct_target']} direct-target-supported drugs. AP is computed on a 1:1 positive/unknown sample, not the full-catalog retrieval task or therapeutic efficacy. Drug popularity is computed on training indications only.

## Demonstrated indication-recovery performance

{md_table(aggregate.reset_index())}

These numbers are the measured results; a changed subtype ranking cannot compensate for poor held-out performance. Initialization ranges on one fixed split are not confidence intervals or external validation. Disease-macro AP reduces domination by heavily annotated diseases but is still a sampled annotation-recovery metric. Historical notebook values (GNN AP 0.819, popularity AP 0.894) are context only; reconstructed inputs/software can change exact values.

## Why the original subtype rankings can collapse

The original code ignores `mutation_fraction` and averages broad, overlapping observed-gene unions. It also collapses disease-gene associations and subtype mutation observations into the same aggregation rule. Gene averaging, PPI smoothing, a drug-driven decoder and annotation popularity can all reduce subtype dependence; the experiments below measure their consequences rather than assuming a cause. Both the original and target-only drug pooling still omit action direction and variant function.

Profile construction is prespecified, with no disease-test or desired-drug tuning: (1) every observed mapped gene equally weighted; (2) all positive genes weighted by subtype observation frequency; (3) retain genes observed in at least two patients and at least 10% of the subtype, cap at 200 highest frequencies, and normalize frequency weights. This is descriptive selection, **not burden-adjusted enrichment**. A fixed node-index tie rule is used; cutoff ties are audited. Six-patient GI.HM-SNV remains uncertain. `profile_summary.csv` reports actual gene counts and effective weighted sizes.

## Subtype sensitivity: differentiation is not validation

The table averages the six subtype pairs and available initialization seeds for the target-only model. Spearman uses tied scores; deterministic node-index tie breaking is used for top20 membership, so rank-list differences near ties can be arbitrary.

{md_table(sensitivity.reset_index())}

For the following table, each condition is compared with that same subtype's selective-frequency profile using identical model weights:

{md_table(changes.reset_index())}

Controls: identical pooled profiles must produce identical scores; zero profiles remove subtype molecular input; a fixed gene permutation preserves each profile's weights/size but destroys gene identity correspondence; swapping the four profiles must swap predictions exactly. Assertions verify identical-profile and swap behavior. The retrained drug-only condition removes all disease molecular context, testing whether labeled disease gene neighborhoods improve the supervised task. The no-PPI condition is a single-seed graph-structure diagnostic. These are complementary tests: perturbation sensitivity shows computational dependence, whereas predictive improvement over controls would support task-relevant information. Random profiles can also generate different rankings; differentiation alone is not success.

`embedding_diagnostics.csv` compares subtype embedding cosines and between-subtype score RMS relative to across-drug score variation. `scores_*.npz` retain full candidate scores for auditing. `hypothesis_rankings.csv` provides raw logits and a separately labeled within-drug subtype contrast (score minus that drug's mean over four subtypes). Contrast is a diagnostic to remove a shared drug preference, **not a second validated predictor**. Hypotheses are not calibrated probabilities, clinical recommendations, or evidence of differential response. No independent subtype response labels or prospective novelty validation are available.

## Interpretation and next steps

Use the measured performance and controls together. If weighting/selectivity changes rankings but not biological predictive evidence, report improved controllability/sensitivity only. If drug-only or popularity matches/exceeds the GNN, acknowledge that the present experiment has not demonstrated a GNN advantage. If permuted profiles produce comparable differentiation, do not infer that differentiated rankings are biologically meaningful. A no-PPI comparison at one seed is exploratory and cannot establish a reliable graph benefit.

The defensible proof of concept is traceable molecular input -> graph model -> explicitly unvalidated association hypotheses, with honest baselines and sensitivity audits. Before selecting hypotheses for further investigation, inspect variant consequences, driver/passenger context, drug action, known CRC indications, assay evidence and source independence. Disease-level evaluation does not validate the four unlabeled subtype predictions. More distinct lists are not a scientific objective on their own.

## Reproduction and preservation

`python3 analysis/crc_experiments.py` creates a new timestamped directory; `--output` must name a nonexistent directory. Defaults are 40 epochs and model seeds 42,43,44 on CPU with two threads and deterministic algorithms. The script never executes notebook cells, downloads data itself, changes original outputs, or writes model results into the original prediction directories. All protected inputs, notebooks and preexisting report/artifact files were fingerprinted before and after and remained unchanged. Versions, hashes, split IDs, sampled pairs, training curves and diagnostic scores are saved. Dependencies are existing pandas, NumPy, SciPy, scikit-learn and PyTorch.
"""
    (out/"crc_experiment_report.md").write_text(report)
    from crc_report import write_summary
    write_summary(out)


if __name__=="__main__":
    main()
