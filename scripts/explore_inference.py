"""Reproduce coverage analysis and masked-edge ranking baselines (no ML).

Run after scripts/prepare_primekg.py. Outputs small reports, not a second graph.
"""

import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT = PROJECT_ROOT / "reports/inference"
SEEDS = (17, 29, 43)
SELECT = {
    "phenotype_protein", "disease_protein", "disease_phenotype_positive",
    "disease_phenotype_negative", "protein_protein", "pathway_protein",
}


def adjacency(pairs):
    result = defaultdict(set)
    for a, b in pairs:
        result[int(a)].add(int(b))
    return result


def expected_metrics(scores, target, excluded, harmonic):
    """Expected metrics over uniform random permutations within score ties."""
    allowed = np.ones(len(scores), dtype=bool)
    allowed[list(excluded - {target})] = False
    values = scores[allowed]
    value = scores[target]
    higher = int(np.count_nonzero(values > value))
    tied = int(np.count_nonzero(values == value))
    assert tied > 0
    rr = (harmonic[higher + tied] - harmonic[higher]) / tied
    return [rr, *[max(0, min(tied, k - higher)) / tied for k in (1, 10, 100)]]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    nodes = pd.read_csv(PROJECT_ROOT / "data/processed/primekg/nodes.csv", dtype={"id": "string"},
                        keep_default_na=False).set_index("index")
    # Canonicalize exact external-ID aliases before splitting and scoring.
    canonical = nodes.reset_index().groupby(["type", "source", "id"])["index"].transform("min")
    mapping = dict(zip(nodes.index, canonical.astype(int)))
    representative = nodes.loc[sorted(set(mapping.values()))]
    types = nodes.type.to_dict()
    pieces, inventory = [], Counter()
    columns = ["x_index", "y_index", "relation", "display_relation", "x_type", "y_type"]
    for chunk in pd.read_csv(PROJECT_ROOT / "data/raw/primekg/kg.csv", usecols=columns, keep_default_na=False, chunksize=250_000):
        inventory.update(chunk.groupby(["relation", "display_relation", "x_type", "y_type"]).size().to_dict())
        pieces.append(chunk[chunk.relation.isin(SELECT)].copy())
    raw = pd.concat(pieces, ignore_index=True)
    old = pd.read_csv(PROJECT_ROOT / "reports/primekg/relationships.csv", keep_default_na=False)
    assert inventory == Counter({tuple(r[:4]): r[4] for r in old.itertuples(index=False, name=None)})
    assert raw.x_type.eq(raw.x_index.map(types)).all()
    assert raw.y_type.eq(raw.y_index.map(types)).all()

    def pairs(relation, first, second):
        subset = raw[raw.relation.eq(relation)]
        forward = subset[subset.x_type.eq(first) & subset.y_type.eq(second)]
        reverse = subset[subset.x_type.eq(second) & subset.y_type.eq(first)]
        a = set(zip(forward.x_index, forward.y_index))
        b = set(zip(reverse.y_index, reverse.x_index))
        if first != second:
            assert a == b, f"Missing reciprocal rows: {relation}"
        else:
            assert a == {(y, x) for x, y in a}
        return {(mapping[x], mapping[y]) for x, y in a}

    gp = pairs("phenotype_protein", "gene/protein", "effect/phenotype")
    gd = pairs("disease_protein", "gene/protein", "disease")
    dp = pairs("disease_phenotype_positive", "disease", "effect/phenotype")
    dn = pairs("disease_phenotype_negative", "disease", "effect/phenotype")
    ppi = pairs("protein_protein", "gene/protein", "gene/protein")
    pathway = pairs("pathway_protein", "gene/protein", "pathway")
    processed = pd.read_csv(PROJECT_ROOT / "data/processed/primekg/phenotype_disease_gene_edges.csv")
    prior = processed[processed.relation.eq("phenotype_protein")]
    assert gp == {(mapping[x], mapping[y]) if types[x] == "gene/protein" else
                  (mapping[y], mapping[x]) for x, y in zip(prior.x_index, prior.y_index)}
    g_to_d, d_to_p = adjacency(gd), adjacency(dp)
    g_to_p, neighbors = adjacency(gp), adjacency(ppi)
    g_to_path, path_to_g = adjacency(pathway), adjacency((p, g) for g, p in pathway)
    genes = sorted(representative[representative.type.eq("gene/protein")].index)
    phenotypes = sorted(representative[representative.type.eq("effect/phenotype")].index)
    ppos = {p: i for i, p in enumerate(phenotypes)}
    indirect = {g: set().union(*(d_to_p[d] for d in g_to_d[g])) for g in genes}
    direct_set = set(g_to_p)
    indirect_set = {g for g in genes if indirect[g]}
    coverage = pd.DataFrame([{
        "gene_index": g, "name": nodes.loc[g, "name"], "external_id": nodes.loc[g, "id"],
        "direct_phenotypes": len(g_to_p[g]), "indirect_positive_phenotypes": len(indirect[g]),
        "diseases": len(g_to_d[g]), "ppi_neighbors": len(neighbors[g]),
        "pathways": len(g_to_path[g]),
    } for g in genes])
    coverage.to_csv(OUT / "gene_coverage.csv", index=False)
    direct_reachable = sum(p in indirect[g] for g, p in gp)
    summary = {
        "canonical_genes": len(genes), "canonical_phenotypes": len(phenotypes),
        "canonical_diseases": int(representative.type.eq("disease").sum()),
        "direct_pairs": len(gp), "direct_genes": len(direct_set),
        "direct_phenotypes": len({p for g, p in gp}),
        "direct_gene_degree_quantiles": coverage.loc[coverage.direct_phenotypes.gt(0), "direct_phenotypes"].quantile([0,.25,.5,.75,1]).to_dict(),
        "indirect_pairs": sum(map(len, indirect.values())),
        "indirect_genes": len(indirect_set), "direct_or_indirect_genes": len(direct_set | indirect_set),
        "indirect_only_genes": len(indirect_set - direct_set),
        "direct_only_genes": len(direct_set - indirect_set),
        "direct_pairs_reachable_via_disease": direct_reachable,
        "gd_pairs": len(gd), "dp_positive_pairs": len(dp), "dp_negative_pairs": len(dn),
        "dp_positive_negative_overlap": len(dp & dn),
        "ppi_unordered_pairs": len({tuple(sorted((g,h))) for g,h in ppi if g != h}),
        "ppi_canonical_self_loops": sum(g == h for g,h in ppi),
        "gene_pathway_pairs": len(pathway), "genes_with_pathways": len({g for g,p in pathway}),
        "genes_with_diseases": len({g for g,d in gd}),
        "genes_with_ppi": len({g for g,h in ppi}),
    }
    examples = []
    for label, candidates in (("documented_and_path", sorted((g,p) for g,p in gp if p in indirect[g])),
                               ("path_only", ((g,p) for g in genes for p in sorted(indirect[g] - g_to_p[g])))):
        for g,p in candidates:
            ds = sorted(d for d in g_to_d[g] if p in d_to_p[d])
            if ds:
                d = ds[0]
                examples.append({"kind": label, "gene": int(g), "gene_name": nodes.loc[g,"name"],
                    "disease": int(d), "disease_name": nodes.loc[d,"name"],
                    "phenotype": int(p), "phenotype_name": nodes.loc[p,"name"]})
            if sum(e["kind"] == label for e in examples) >= 3:
                break
    (OUT / "examples.json").write_text(json.dumps(examples, indent=2) + "\n")

    # Sparse score accumulation; only the final candidate score vector is dense.
    def vector(counter):
        result = np.zeros(len(phenotypes), dtype=float)
        for p, value in counter.items():
            result[ppos[p]] = value
        return result

    disease_scores = {}
    for g in direct_set:
        counts = Counter()
        for d in g_to_d[g]:
            # Penalize diseases with many annotations, without assuming frequency.
            weight = 1 / len(d_to_p[d]) if d_to_p[d] else 0
            for p in d_to_p[d]:
                counts[p] += weight
        disease_scores[g] = vector(counts)
    harmonic = np.r_[0., np.cumsum(1 / np.arange(1, len(phenotypes) + 1))]
    pair_list = sorted(gp)
    rows, splits, audits = [], [], []
    for mode in ("edge_holdout", "gene_annotation_holdout"):
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            if mode == "edge_holdout":
                chosen = rng.choice(len(pair_list), size=round(.2 * len(pair_list)), replace=False)
                test = {pair_list[i] for i in chosen}
            else:
                annotated = sorted(direct_set)
                chosen = rng.choice(len(annotated), size=round(.2 * len(annotated)), replace=False)
                held_genes = {annotated[i] for i in chosen}
                test = {(g,p) for g,p in gp if g in held_genes}
            train = gp - test
            assert not train & test
            train_gp = adjacency(train)
            if mode == "gene_annotation_holdout":
                assert not set(train_gp) & held_genes
            for g,p in sorted(test):
                splits.append({"mode": mode, "seed": seed, "gene_index": g, "phenotype_index": p})
            popularity = vector(Counter(p for g,p in train))
            path_counts = {}
            for path, members in path_to_g.items():
                path_counts[path] = Counter(p for g in members for p in train_gp[g])
            audits.append({"mode": mode, "seed": seed, "train_pairs": len(train),
                "test_pairs": len(test), "test_genes": len({g for g,p in test}),
                "test_phenotypes_without_training_labels": len({p for g,p in test} - {p for g,p in train}),
                "test_pairs_with_disease_path": sum(p in indirect[g] for g,p in test)})
            for g, targets in sorted(adjacency(test).items()):
                ppi_count = Counter(p for h in neighbors[g] if h != g for p in train_gp[h])
                pathway_count = Counter()
                for path in g_to_path[g]:
                    # Exclude the query gene's own annotations from shared-pathway votes.
                    weight = 1 / max(1, len(path_to_g[path]) - 1)
                    for p, count in path_counts[path].items():
                        pathway_count[p] += weight * (count - int(p in train_gp[g]))
                ppi_vec, path_vec = vector(ppi_count), vector(pathway_count)
                scores = {"direct_lookup": vector(Counter(train_gp[g])),
                    "phenotype_popularity": popularity,
                    "disease_paths": disease_scores[g],
                    "ppi_transfer": ppi_vec, "shared_pathway_transfer": path_vec}
                # Fixed equal channel weighting, not tuned on held-out labels.
                ppi_norm = ppi_vec / max(1., ppi_vec.max())
                path_norm = path_vec / max(1., path_vec.max())
                scores["ppi_pathway_equal_mix"] = (ppi_norm + path_norm) / 2
                for method, score in scores.items():
                    excluded = {ppos[p] for p in g_to_p[g]}
                    for p in sorted(targets):
                        values = expected_metrics(score, ppos[p], excluded, harmonic)
                        rows.append({"mode": mode, "seed": seed, "gene_index": g,
                            "phenotype_index": p, "method": method,
                            **dict(zip(["mrr", "hits1", "hits10", "hits100"], values)),
                            "nonzero_target": bool(score[ppos[p]] > 0),
                            "disease_path_available": p in indirect[g]})
    results = pd.DataFrame(rows)
    keys = ["mode", "seed", "method"]
    metrics = ["mrr", "hits1", "hits10", "hits100", "nonzero_target"]
    aggregate = results.groupby(keys)[metrics].mean().reset_index()
    macro = results.groupby(keys + ["gene_index"]).mrr.mean().groupby(keys).mean().rename("gene_macro_mrr")
    aggregate = aggregate.merge(macro.reset_index(), on=keys)
    aggregate.to_csv(OUT / "baseline_metrics.csv", index=False)
    results.groupby(["mode", "method", "disease_path_available"])[metrics].mean().reset_index().to_csv(
        OUT / "baseline_path_strata.csv", index=False)
    pd.DataFrame(splits).to_csv(OUT / "test_splits.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "split_audit.csv", index=False)
    digest = hashlib.sha256()
    with open(PROJECT_ROOT / "data/raw/primekg/kg.csv", "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    summary.update({"raw_sha256": digest.hexdigest(), "seeds": SEEDS,
        "pandas_version": pd.__version__, "numpy_version": np.__version__})
    (OUT / "coverage_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(aggregate.groupby(["mode", "method"])[metrics + ["gene_macro_mrr"]].mean().to_string())


if __name__ == "__main__":
    main()
