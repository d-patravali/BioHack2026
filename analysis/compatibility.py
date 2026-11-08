"""Local ECoGraph/PrimeKG compatibility analysis; never merges source graphs."""

import argparse
import hashlib
import json
import re
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import hypergeom, norm

SUBTYPES = ["GI.CIN", "GI.HM-indel", "GI.HM-SNV", "GI.GS"]
FILES = ["vertices_gene.csv", "vertices_subtype.csv", "vertices_sample.csv", "vertices_age.csv",
         "edges_has_gene_part001.csv", "edges_has_subtype.csv", "edges_in_pathway.csv", "edges_is_of_age.csv"]
NODE_FIELDS = ["index", "id", "type", "name", "source"]


def read_csv(path):
    try:
        return pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pct(numerator, denominator):
    return 100 * numerator / denominator if denominator else np.nan


def table(df):
    """Markdown without an optional tabulate dependency."""
    def fmt(value):
        if pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{value:.2f}"
        return str(value).replace("|", "\\|")
    return "\n".join(["| " + " | ".join(df.columns) + " |",
                      "| " + " | ".join(["---"] * len(df.columns)) + " |"] +
                     ["| " + " | ".join(map(fmt, row)) + " |" for row in df.itertuples(index=False, name=None)])


def bh(pvalues):
    """Benjamini–Hochberg correction, including p=1 untestable hypotheses."""
    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues, kind="stable")
    ranked = pvalues[order] * len(pvalues) / np.arange(1, len(pvalues) + 1)
    corrected = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
    result = np.empty(len(pvalues))
    result[order] = corrected
    return result


def inspect_sources(paths, out):
    manifest, lines = [], []
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"Required source missing: {path}")
        frame = pd.read_csv(path, nrows=3, dtype=str, keep_default_na=False) if path.stat().st_size > 2 else read_csv(path)
        info = {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path),
                "columns": list(frame.columns), "representative_records": frame.to_dict("records")}
        manifest.append(info)
        lines.extend([str(path), f"  bytes={info['bytes']}", f"  columns={info['columns']}",
                      "  records=" + json.dumps(info["representative_records"])])
    (out / "source_inspection.txt").write_text("\n".join(lines) + "\n")
    (out / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def ecograph_audit(frames, out):
    vertices = {name: set(d["~id"]) for name, d in frames.items() if name.startswith("vertices")}
    declared = set().union(*vertices.values())
    rows, dangling = [], []
    for name, frame in frames.items():
        row = {"file": name, "rows": len(frame), "columns": "|".join(frame.columns),
               "duplicate_rows": int(frame.duplicated().sum()),
               "blank_cells": int(frame.eq("").sum().sum()),
               "missing_marker_cells": int(frame.isin(["NA", "N/A", "NaN", "null"]).sum().sum()),
               "duplicate_ids": int(frame["~id"].duplicated().sum()) if "~id" in frame else 0,
               "blank_ids": int(frame["~id"].eq("").sum()) if "~id" in frame else 0,
               "labels": "|".join(sorted(frame["~label"].unique())) if "~label" in frame else "unavailable",
               "status": "empty_no_schema" if len(frame.columns) == 0 else "parsed"}
        if "~from" in frame:
            expected = {
                "edges_has_subtype.csv": ("vertices_sample.csv", "vertices_subtype.csv", "HAS_SUBTYPE"),
                "edges_is_of_age.csv": ("vertices_sample.csv", "vertices_age.csv", "IS_OF_AGE"),
                "edges_in_pathway.csv": ("vertices_gene.csv", None, "IN_PATHWAY"),
                "edges_has_gene_part001.csv": ("vertices_sample.csv", "vertices_gene.csv", "HAS_GENE"),
            }[name]
            wrong = ~frame["~from"].isin(vertices[expected[0]]) | frame["~label"].ne(expected[2])
            wrong |= ~frame["~to"].isin(vertices[expected[1]]) if expected[1] else ~frame["~to"].str.startswith("pathway:")
            row["endpoint_or_label_type_mismatch_rows"] = int(wrong.sum())
            if wrong.any():
                raise ValueError(f"Unexpected endpoint types/direction/label in {name}")
            for col in ("~from", "~to"):
                values = frame.loc[~frame[col].isin(declared), col]
                row[f"undeclared_{col[1:]}_edge_rows"] = len(values)
                row[f"undeclared_{col[1:]}_unique_nodes"] = values.nunique()
                dangling.extend({"file": name, "endpoint": col, "id": v} for v in sorted(values.unique()))
        rows.append(row)
    pd.DataFrame(rows).to_csv(out / "ecograph_integrity.csv", index=False)
    pd.DataFrame(dangling, columns=["file", "endpoint", "id"]).to_csv(out / "undeclared_endpoints.csv", index=False)
    return pd.DataFrame(rows)


def scan_primekg(path, out):
    """Use the existing preprocessing schema; rebuild metadata from raw chunks."""
    parts, relation_counts, edge_parts = [], Counter(), []
    blanks = Counter()
    raw_rows = 0
    # All metadata columns are needed for identifier/alias auditing.
    columns = ["relation", "display_relation"] + [f"{s}_{f}" for s in ("x", "y") for f in NODE_FIELDS]
    for chunk in pd.read_csv(path, usecols=columns, dtype=str, keep_default_na=False, chunksize=250_000):
        raw_rows += len(chunk)
        blanks.update(chunk.eq("").sum().to_dict())
        relation_counts.update(chunk.groupby(["relation", "display_relation", "x_type", "y_type"]).size().to_dict())
        for side in ("x", "y"):
            nodes = chunk[[f"{side}_{field}" for field in NODE_FIELDS]].copy()
            nodes.columns = NODE_FIELDS
            parts.append(nodes.drop_duplicates())
        # Preserve orientation and both labels; retain only relevant drug edges.
        selected = chunk[chunk.relation.isin(["drug_protein", "indication", "contraindication", "off-label use"])]
        edge_parts.append(selected[["x_index", "y_index", "relation", "display_relation", "x_type", "y_type"]].copy())
    observations = pd.concat(parts, ignore_index=True).drop_duplicates()
    conflicts = observations[observations.duplicated("index", keep=False)]
    conflicts.to_csv(out / "primekg_metadata_conflicts.csv", index=False)
    if len(conflicts):
        raise ValueError("PrimeKG index metadata conflicts require resolution")
    aliases = observations[observations.duplicated(["type", "source", "id"], keep=False)]
    aliases.sort_values(["type", "source", "id", "index"]).to_csv(out / "primekg_duplicate_identifiers.csv", index=False)
    nodes = observations.sort_values("index").reset_index(drop=True)
    assert nodes["index"].is_unique
    edges = pd.concat(edge_parts, ignore_index=True).drop_duplicates()
    assert set(edges.x_index) | set(edges.y_index) <= set(nodes["index"])
    inventory = pd.DataFrame([(*key, value) for key,value in sorted(relation_counts.items())],
        columns=["relation", "display_relation", "x_type", "y_type", "raw_rows"])
    inventory.to_csv(out / "primekg_relations.csv", index=False)
    nodes.groupby(["type", "source"]).size().rename("nodes").reset_index().to_csv(out / "primekg_node_types.csv", index=False)
    (out / "primekg_summary.json").write_text(json.dumps({"raw_rows": raw_rows, "node_indices": len(nodes),
        "relation_labels": inventory.relation.nunique(), "blank_cells": dict(blanks),
        "repeated_external_identifier_rows": len(aliases)}, indent=2) + "\n")
    return nodes, edges, inventory


def map_genes(eco, prime):
    """Unique exact symbol lookup to stable IDs; only whitespace is normalized."""
    by_name = prime.groupby("name")
    result = []
    for row in eco.itertuples(index=False, name=None):
        vertex, _, symbol = row[:3]
        invalid = symbol == "patient_id"
        key = symbol if symbol in by_name.groups else symbol.strip()
        candidates = by_name.get_group(key) if key in by_name.groups else prime.iloc[:0]
        ids = sorted(set(candidates.id))
        status = "excluded_metadata_artifact" if invalid else (
            "ambiguous" if len(ids) > 1 else "exact_symbol" if len(ids) == 1 and key == symbol else
            "normalized_whitespace_symbol" if len(ids) == 1 else "unmatched")
        result.append({"entity_type": "gene", "ecograph_id": vertex, "ecograph_label": symbol,
            "status": status, "candidate_primekg_ids": "|".join(ids) if not invalid else "",
            "primekg_id": ids[0] if len(ids) == 1 and not invalid else "",
            "primekg_source": "NCBI" if len(ids) == 1 and not invalid else "",
            "primekg_indices": "|".join(sorted(candidates["index"])) if not invalid else "",
            "mapping_basis": "exact unique symbol to existing NCBI ID; not independently verified stable-ID crosswalk"})
    return pd.DataFrame(result)


def pathway_audit(path_edges, prime, out):
    def normalize(text):
        return re.sub(r"[^A-Z0-9]+", "_", text.upper()).strip("_")
    lookup = {}
    for row in prime.itertuples(index=False):
        lookup.setdefault(normalize(row.name), set()).add(row.id)
    rows = []
    for value in sorted(path_edges["~to"].unique()):
        label = value.removeprefix("pathway:")
        collection = label.split("_", 1)[0]
        # Names are candidate suggestions only, even when normalized labels coincide.
        candidates = lookup.get(normalize(label.removeprefix("REACTOME_")), set()) if collection == "REACTOME" else set()
        stable = prime[prime.id.eq(label)]
        confident = sorted(set(stable.id))
        rows.append({"entity_type": "pathway", "ecograph_id": value, "ecograph_label": label,
            "collection": collection, "status": "exact_stable_id" if len(confident) == 1 else
            "ambiguous_label_candidate" if len(candidates) > 1 else
            "unverified_label_candidate" if candidates else "unmatched_no_crosswalk",
            "candidate_primekg_ids": "|".join(sorted(candidates)),
            "primekg_id": confident[0] if len(confident) == 1 else "",
            "primekg_source": "REACTOME" if confident else "",
            "mapping_basis": "stable-ID equality only; collection-name normalization is not accepted equivalence"})
    audit = pd.DataFrame(rows)
    audit.to_csv(out / "pathway_mapping_audit.csv", index=False)
    summary = pd.DataFrame([{"ecograph_pathways": len(audit), "primekg_pathways": prime.id.nunique(),
        "confident_matches": audit.primekg_id.ne("").sum(),
        "confident_mapping_pct": pct(audit.primekg_id.ne("").sum(), len(audit)),
        "unmatched_without_confident_mapping": audit.primekg_id.eq("").sum(),
        "ambiguous_label_candidates": audit.status.eq("ambiguous_label_candidate").sum(),
        "unique_unverified_label_candidates": audit.status.eq("unverified_label_candidate").sum(),
        "biological_overlap_pct": np.nan,
        "interpretation": "No stable pathway IDs/crosswalk supplied; confident mapping coverage is not biological overlap"}])
    summary.to_csv(out / "pathway_overlap_summary.csv", index=False)
    return audit, summary


def therapeutic_context(edges, nodes, out):
    metadata = nodes.set_index("index")
    drug_gene = edges[edges.relation.eq("drug_protein")]
    forward = drug_gene[drug_gene.x_type.eq("drug")]
    reverse = drug_gene[drug_gene.y_type.eq("drug")]
    fields = ["x_index", "y_index", "display_relation"]
    assert set(map(tuple, forward[fields].values)) == set(zip(reverse.y_index, reverse.x_index, reverse.display_relation))
    drug_gene = forward.assign(drug_id=forward.x_index.map(metadata.id),
        gene_id=forward.y_index.map(metadata.id)).drop_duplicates(["drug_id", "gene_id", "display_relation"])
    drug_gene.groupby("display_relation").size().rename("unique_drug_gene_pairs").reset_index().to_csv(out / "drug_gene_relation_counts.csv", index=False)
    targets = drug_gene[drug_gene.display_relation.eq("target")]
    targets[["drug_id", "gene_id", "display_relation"]].to_csv(out / "primekg_direct_targets.csv", index=False)
    indication = edges[edges.relation.eq("indication") & edges.x_type.eq("drug")].copy()
    reverse_indication = edges[edges.relation.eq("indication") & edges.y_type.eq("drug")]
    assert set(zip(indication.x_index, indication.y_index)) == set(zip(reverse_indication.y_index, reverse_indication.x_index))
    indication["drug_id"] = indication.x_index.map(metadata.id)
    indication["disease_id"] = indication.y_index.map(metadata.id)
    indication = indication.drop_duplicates(["drug_id", "disease_id"])
    disease = nodes[nodes.type.eq("disease")]
    mask = disease.name.str.contains(r"colorectal|\bcolon\b|\brectal\b", case=False, regex=True)
    crc = disease[mask].copy()
    counts = indication.groupby("y_index").drug_id.nunique()
    crc["known_indication_drugs"] = crc["index"].map(counts).fillna(0).astype(int)
    crc["search_note"] = "lexical candidate, not an asserted mapping to ECoGraph genomic subtype"
    crc.to_csv(out / "primekg_colorectal_disease_candidates.csv", index=False)
    crc_indications = indication[indication.y_index.isin(crc["index"])].copy()
    crc_indications["disease_name"] = crc_indications.y_index.map(metadata.name)
    crc_indications["drug_name"] = crc_indications.x_index.map(metadata.name)
    crc_indications[["drug_id", "drug_name", "disease_id", "disease_name"]].to_csv(out / "primekg_crc_indications.csv", index=False)
    info = {"direct_target_pairs": len(targets), "target_genes": targets.gene_id.nunique(),
        "target_drugs": targets.drug_id.nunique(), "indication_pairs": len(indication),
        "indication_drugs": indication.drug_id.nunique(), "indication_diseases": indication.disease_id.nunique(),
        "crc_lexical_candidates": len(crc), "crc_indication_pairs": len(crc_indications),
        "crc_indication_drugs": crc_indications.drug_id.nunique()}
    return targets, info, crc


def supplementary_matrix(base, frames, genes, out):
    """Explicitly separate recovery analysis; never writes source exports."""
    path = base.parent / "coad_mut_data.csv"
    subtype_path = base.parent / "coad_subtype_data.csv"
    if not path.is_file() or not subtype_path.is_file():
        return None, None, {"status": "unavailable"}
    matrix = pd.read_csv(path, keep_default_na=False)
    assert {"patient_id", "Tumor_Sample_Barcode"} <= set(matrix.columns)
    assert matrix.patient_id.is_unique
    assert matrix.patient_id.eq(matrix.Tumor_Sample_Barcode.str[:12]).all()
    sample_ids = "sample:" + matrix.patient_id
    assert set(sample_ids) == set(frames["vertices_sample.csv"]["~id"])
    gene_cols = [c for c in matrix.columns if c not in {"patient_id", "Tumor_Sample_Barcode"}]
    assert set(gene_cols) == set(genes.loc[genes.status.ne("excluded_metadata_artifact"), "ecograph_label"])
    x = matrix[gene_cols].to_numpy()
    assert np.isin(x, [0, 1]).all(), "Recovery requires observed binary matrix, no implicit coercion"
    x = pd.DataFrame(x.astype(np.uint8), index=sample_ids, columns=gene_cols)
    sub_edges = frames["edges_has_subtype.csv"].drop_duplicates(["~from", "~to"])
    assert sub_edges["~from"].is_unique, "Multi-subtype assignments require explicit policy"
    annotation = sub_edges.set_index("~from")["~to"].str.removeprefix("subtype:")
    source_sub = read_csv(subtype_path)
    assert set(zip(sub_edges["~from"], sub_edges["~to"])) <= set(zip(
        "sample:" + source_sub["pan.samplesID"], "subtype:" + source_sub.Subtype_Selected))
    assert set(annotation) == set(SUBTYPES)
    burden = x.sum(axis=1).rename("observed_mutated_genes")
    burden_frame = burden.to_frame().join(annotation.rename("subtype"))
    burden_frame.index.name = "ecograph_sample_id"
    burden_frame.to_csv(out / "supplementary_sample_burden.csv")
    info = {"status": "verified_supplementary_reconstruction", "matrix_path": str(path),
        "samples": len(x), "genes": len(gene_cols), "observations_all_samples": int(x.values.sum()),
        "typed_samples": len(annotation), "untyped_samples": len(x) - len(annotation),
        "matrix_sha256": sha256(path), "subtype_sha256": sha256(subtype_path),
        "duplicate_patient_records": int(matrix.patient_id.duplicated().sum())}
    return x, annotation, info


def exploratory_enrichment(x, labels, mapping, targets, out):
    """Descriptive screens, not validated subtype-specific biological markers.

    Fisher one-sided tests compare prevalence; pooled burden quintiles give a
    separate one-sided CMH normal approximation, with BH across all subtype/gene
    hypotheses for each test family. Tiny strata/groups preclude strong claims.
    """
    x = x.loc[labels.index]
    burden = x.sum(axis=1)
    strata = pd.qcut(burden, q=5, duplicates="drop")
    output = []
    for subtype in SUBTYPES:
        member = labels.eq(subtype).to_numpy()
        n1, n0 = int(member.sum()), int((~member).sum())
        a = x.loc[member].sum().to_numpy(dtype=float)
        c = x.loc[~member].sum().to_numpy(dtype=float)
        fisher = hypergeom.sf(a - 1, n1+n0, a+c, n1)
        delta, variance = np.zeros(x.shape[1]), np.zeros(x.shape[1])
        supported1 = supported0 = mixed_strata = 0
        for bucket in strata.cat.categories:
            selected = strata.eq(bucket).to_numpy()
            one, zero = selected & member, selected & ~member
            sn1, sn0 = int(one.sum()), int(zero.sum())
            if sn1 == 0 or sn0 == 0:
                continue
            mixed_strata += 1
            supported1 += sn1
            supported0 += sn0
            sa = x.loc[one].sum().to_numpy(dtype=float)
            sm = x.loc[selected].sum().to_numpy(dtype=float)
            sn = sn1 + sn0
            delta += sa - sn1 * sm / sn
            variance += sn1 * sn0 * sm * (sn-sm) / (sn*sn*(sn-1))
        z = np.divide(delta, np.sqrt(variance), out=np.zeros_like(delta), where=variance > 0)
        cmh = np.where(variance > 0, norm.sf(z), 1.)
        # Minimal screening guardrails, not a guarantee of asymptotic validity.
        usable = (n1 >= 20 and n0 >= 20 and supported1 / n1 >= .5 and supported0 / n0 >= .5)
        output.append(pd.DataFrame({"subtype": subtype, "symbol": x.columns,
            "subtype_samples": n1, "other_typed_samples": n0, "subtype_observed": a.astype(int),
            "other_observed": c.astype(int), "subtype_frequency": a/n1, "other_frequency": c/n0,
            "fisher_one_sided_p": fisher, "burden_cmh_normal_p_unfiltered": cmh,
            "burden_cmh_one_sided_p": cmh if usable else np.ones_like(cmh),
            "burden_screen_eligible": usable,
            "mixed_burden_strata": mixed_strata, "subtype_samples_in_mixed_strata": supported1,
            "other_samples_in_mixed_strata": supported0,
            "interpretation": "exploratory association screen; not validated enrichment"}))
    result = pd.concat(output, ignore_index=True)
    result["fisher_bh_q"] = bh(result.fisher_one_sided_p)
    result["burden_cmh_bh_q"] = bh(result.burden_cmh_one_sided_p)
    result["screen_candidate"] = (result.fisher_bh_q.lt(.05) & result.burden_cmh_bh_q.lt(.05)
                                  & result.subtype_frequency.gt(result.other_frequency))
    gene_map = mapping.set_index("ecograph_label").primekg_id.to_dict()
    result["primekg_gene_id"] = result.symbol.map(gene_map).fillna("")
    result["direct_drug_target"] = result.primekg_gene_id.isin(set(targets.gene_id))
    result.to_csv(out / "supplementary_enrichment_screens.csv", index=False)
    rows = []
    for subtype, d in result.groupby("subtype", sort=False):
        screen = d[d.screen_candidate]
        usable = bool(d.burden_screen_eligible.iloc[0])
        rows.append({"subtype": subtype, "raw_fisher_bh_hits": int(d.fisher_bh_q.lt(.05).sum()),
            "burden_screen_eligible": usable,
            "subtype_samples_in_mixed_strata": int(d.subtype_samples_in_mixed_strata.iloc[0]),
            "other_samples_in_mixed_strata": int(d.other_samples_in_mixed_strata.iloc[0]),
            "burden_cmh_bh_hits": int(d.burden_cmh_bh_q.lt(.05).sum()) if usable else np.nan,
            "joint_exploratory_screen_candidates": len(screen) if usable else np.nan,
            "screen_candidates_mapped": int(screen.primekg_gene_id.ne("").sum()) if usable else np.nan,
            "screen_candidates_mapped_pct": pct(screen.primekg_gene_id.ne("").sum(), len(screen)),
            "screen_candidates_with_direct_targets": int(screen.direct_drug_target.sum()) if usable else np.nan,
            "validated_enriched_genes": np.nan,
            "interpretation": "No validated enriched set: asymptotic burden screen, sparse groups, no external replication"})
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "supplementary_enrichment_summary.csv", index=False)
    return result, summary


def compare_sets(sets, basis):
    rows = []
    for a, b in combinations(SUBTYPES, 2):
        union = sets[a] | sets[b]
        rows.append({"basis": basis, "subtype_a": a, "subtype_b": b,
            "genes_a": len(sets[a]), "genes_b": len(sets[b]),
            "overlap": len(sets[a] & sets[b]), "union": len(union),
            "jaccard": len(sets[a] & sets[b]) / len(union) if union else np.nan})
    return pd.DataFrame(rows)


def subtype_analysis(frames, mapping, targets, x, labels, out):
    folder = out / "subtype_gene_sets"
    folder.mkdir(exist_ok=True)
    observed = {}
    gene_rows, target_rows, burden_rows = [], [], []
    for subtype in SUBTYPES:
        n = int(labels.eq(subtype).sum()) if labels is not None else int(frames["edges_has_subtype.csv"]["~to"].eq("subtype:" + subtype).sum())
        counts = x.loc[labels.index[labels.eq(subtype)]].sum() if x is not None else pd.Series(dtype=float)
        counts = counts[counts.gt(0)].astype(int)
        observed[subtype] = set(counts.index)
        data = pd.DataFrame({"symbol": counts.index, "unique_observed_samples": counts.values})
        data["sample_frequency"] = data.unique_observed_samples / n
        data["basis"] = "supplementary_binary_mutation_matrix" if x is not None else "unavailable_empty_sample_gene_export"
        data = data.merge(mapping, left_on="symbol", right_on="ecograph_label", how="left", validate="one_to_one")
        data["direct_drug_target"] = data.primekg_id.isin(set(targets.gene_id))
        data.to_csv(folder / f"{subtype}.csv", index=False)
        matched = data.primekg_id.ne("")
        target_genes = set(data.loc[matched, "primekg_id"]) & set(targets.gene_id)
        target_links = targets[targets.gene_id.isin(target_genes)]
        gene_rows.append({"basis": data.basis.iloc[0] if len(data) else "unavailable", "subtype": subtype,
            "samples": n, "associated_genes": len(data) if x is not None else np.nan,
            "sample_gene_observations": int(counts.sum()) if x is not None else np.nan,
            "mapped_genes": int(matched.sum()) if x is not None else np.nan,
            "mapping_pct": pct(matched.sum(), len(data)),
            "unmatched_genes": int(data.status.eq("unmatched").sum()) if x is not None else np.nan,
            "ambiguous_genes": int(data.status.eq("ambiguous").sum()) if x is not None else np.nan,
            "validated_enriched_genes": np.nan})
        target_rows.append({"basis": gene_rows[-1]["basis"], "subtype": subtype,
            "mapped_genes": gene_rows[-1]["mapped_genes"],
            "mapped_unique_ncbi_ids": data.loc[matched, "primekg_id"].nunique() if x is not None else np.nan,
            "mapped_genes_with_direct_targets": int(data.direct_drug_target.sum()) if x is not None else np.nan,
            "direct_target_pct_of_mapped_genes": pct(data.direct_drug_target.sum(), matched.sum()),
            "unique_connected_drugs": target_links.drug_id.nunique() if x is not None else np.nan,
            "unique_direct_drug_gene_pairs": len(target_links) if x is not None else np.nan,
            "validated_enriched_genes_with_targets": np.nan})
        if x is not None:
            values = x.loc[labels.index[labels.eq(subtype)]].sum(axis=1)
            burden_rows.append({"subtype": subtype, "samples": n, "minimum": values.min(),
                "q25": values.quantile(.25), "median": values.median(), "q75": values.quantile(.75),
                "maximum": values.max(), "mean": values.mean()})
    common = set.intersection(*observed.values()) if observed else set()
    for row in gene_rows:
        s = row["subtype"]
        others = set.union(*(observed[t] for t in SUBTYPES if t != s))
        row["shared_with_other_subtypes"] = len(observed[s] & others) if x is not None else np.nan
        row["exclusive_observed_genes"] = len(observed[s] - others) if x is not None else np.nan
        file = folder / f"{s}.csv"
        detail = read_csv(file)
        detail["observed_in_other_subtypes"] = detail.symbol.isin(others)
        detail["exclusive_observed_in_this_cohort"] = ~detail.symbol.isin(others)
        detail.to_csv(file, index=False)
    gene_summary, target_summary = pd.DataFrame(gene_rows), pd.DataFrame(target_rows)
    gene_summary.to_csv(out / "subtype_gene_coverage.csv", index=False)
    target_summary.to_csv(out / "subtype_drug_target_coverage.csv", index=False)
    overlap = compare_sets(observed, "supplementary_observed_gene_sets")
    if x is None:
        overlap[["genes_a", "genes_b", "overlap", "union", "jaccard"]] = np.nan
    overlap.to_csv(out / "subtype_gene_overlap.csv", index=False)
    pd.DataFrame({"symbol": sorted(common)}).to_csv(out / "genes_shared_all_subtypes.csv", index=False)
    burden_summary = pd.DataFrame(burden_rows)
    burden_summary.to_csv(out / "supplementary_burden_summary.csv", index=False)
    # Export-only conclusions remain missing, not falsely zero, when no observations exist.
    primary = gene_summary.copy()
    if len(frames["edges_has_gene_part001.csv"]) == 0:
        for c in primary.columns:
            if c not in {"subtype", "samples", "basis"}:
                primary[c] = np.nan
        primary["basis"] = "Neptune_export: sample-gene observations unavailable"
    primary.to_csv(out / "neptune_only_subtype_status.csv", index=False)
    return observed, gene_summary, target_summary, overlap, burden_summary


def build_report(out, eco_audit, nodes, mapping, overall, pathway_summary,
                 sub, drug, overlap, burden, enrichment, context, recovery, crc):
    allshared = pd.read_csv(out / "genes_shared_all_subtypes.csv")
    if recovery.get("status") != "verified_supplementary_reconstruction":
        (out / "compatibility_report.md").write_text(
            "# ECoGraph–PrimeKG export-only compatibility\n\n"
            "## Verdict: VIABLE WITH MODIFICATIONS (conditional technical feasibility)\n\n"
            "The sample–gene edge export is empty with no schema. Subtype gene sets, frequencies, "
            "drug-target coverage and cross-subtype molecular overlap cannot be computed from these exports. "
            "Supplementary recovery was disabled or unavailable; NA values are not zero observations. "
            "The export is not currently merge-ready. No graph or model was built.\n\n"
            "## Overall symbol mapping\n\n" + table(overall) + "\n\n"
            "`patient_id` is an excluded metadata artifact. Exact unique symbols are provisional matches "
            "to NCBI IDs; no aliases are silently resolved.\n\n"
            "## Pathways\n\n" + table(pathway_summary) + "\n\n"
            "Pathway names are not stable-ID crosswalks; biological overlap remains unknown. "
            "Pathway endpoints lack supplied vertices.\n\n"
            "## Therapeutic supervision\n\n"
            f"PrimeKG has {context['direct_target_pairs']} direct drug-target pairs and "
            f"{context['indication_pairs']} indication pairs. These do not supply subtype efficacy labels. "
            "Repair sample–gene exports, audit stable IDs/pathways and define an independently validated "
            "Drug–Disease to Drug–Subtype transfer task before integration or training.\n")
        return
    prime_info = json.loads((out / "primekg_summary.json").read_text())
    prime_counts = nodes.groupby("type").size()
    eco_counts = eco_audit.set_index("file").rows
    n_typed = int(sub.samples.sum())
    n_samples = int(eco_counts["vertices_sample.csv"])
    n_eligible = int(overall.eligible_ecograph_genes.iloc[0])
    dup_nodes = read_csv(out / "primekg_duplicate_identifiers.csv")
    path_audit = read_csv(out / "pathway_mapping_audit.csv")
    collection_counts = path_audit.collection.value_counts().to_dict()
    sample_distribution = ", ".join(f"{int(r.samples)} {r.subtype}" for r in sub.itertuples())
    examples = []
    for symbol in ["APC", "TP53", "KRAS", "BRAF", "PIK3CA"]:
        match = mapping[mapping.ecograph_label.eq(symbol)]
        if not len(match):
            continue
        row = {"symbol": symbol, "NCBI ID": match.primekg_id.iloc[0]}
        target_flag = False
        for subtype in SUBTYPES:
            detail = read_csv(out / "subtype_gene_sets" / f"{subtype}.csv")
            selected = detail[detail.symbol.eq(symbol)]
            row[subtype + " observed %"] = float(selected.sample_frequency.iloc[0]) * 100 if len(selected) else 0.
            target_flag |= bool(len(selected) and selected.direct_drug_target.iloc[0] == "True")
        row["PrimeKG direct target"] = target_flag
        examples.append(row)
    example_frame = pd.DataFrame(examples)
    example_frame.to_csv(out / "selected_gene_examples.csv", index=False)
    report = f"""# ECoGraph–PrimeKG compatibility assessment

## Verdict: VIABLE WITH MODIFICATIONS

**Technical gene-level integration is plausible, but the supplied Neptune export is not ready to merge.** Its sample–gene file is a two-byte blank line with no header or observations. Subtype molecular networks cannot be computed from those eight exports alone. All subtype gene-set numbers below come from a separately verified, local binary mutation matrix, not from the empty edge file. No source files were repaired, no graph was merged, and no GNN, agent, drug ranking, or efficacy prediction was implemented.

The conditional verdict means: repair/validate the sample–gene export, remove a metadata-as-gene artifact, preserve unharmonized pathways separately, and define a defensible prediction/evaluation target. **Predictive feasibility and therapeutic efficacy remain unestablished.**

## Research question and inputs

Can GI.CIN, GI.HM-indel, GI.HM-SNV, and GI.GS molecular observations be linked to PrimeKG for downstream drug–subtype hypothesis generation? This assessment separates technical compatibility, molecular associations, target connectivity, therapeutic supervision, and efficacy. Inputs are `kg.csv` and `ECoGraph/data/neptune_csv/` relative to the project root. `source_inspection.txt` records headers and example records; `source_manifest.json` records sizes and SHA-256 digests.

### ECoGraph export structure and integrity

{table(eco_audit[["file", "rows", "status", "duplicate_ids", "duplicate_rows"]])}

Declared node labels are Gene, Sample, Subtype, and Age. Gene vertices contain `~id`, `~label`, `symbol:String`, `exp_hr_coeff:String`, and `hr_pval:String`. Symbols, not Entrez/Ensembl IDs, define `gene:SYMBOL` identifiers. Sample vertices use `sample:TCGA-...` patient-level IDs plus a patient barcode and sex. Subtypes are `subtype:GI.*`; age IDs are `age:<=50` and `age:>50`. The two survival fields contain literal missing markers, counted separately from blanks in `ecograph_integrity.csv`; they are not subtype evidence or drug-response labels.

Edges use `~id`, `~from`, `~to`, and `~label`. Verified directions are Sample → Subtype (`HAS_SUBTYPE`), Gene → Pathway (`IN_PATHWAY`), and Sample → Age (`IS_OF_AGE`). The sample–gene file has no schema, so its direction and semantics cannot be verified from the export. The accompanying notebook intends Sample → Gene (`HAS_GENE`) with nonzero mutation observations, but this is code intent, not populated edge evidence.

The {n_typed} subtype edges give one assignment each to {sample_distribution} samples. There are {n_samples} sample vertices; {n_samples - n_typed} have no exported subtype. All existing subtype and age edge endpoints resolve to supplied vertices. All pathway source genes resolve; **every pathway destination is undeclared** because no `vertices_pathway.csv` is supplied: {int(eco_counts['edges_in_pathway.csv']):,} edges reference {len(path_audit):,} distinct pathway endpoints. None of the nonempty ECoGraph files has duplicated IDs or rows or blank cells. Missing marker tokens still occur in gene properties.

**Metadata artifact:** gene vertices contain `gene:patient_id`. The mutation matrix has a metadata column `patient_id`, and the notebook's `gene_columns` exclusion list omits it. This is not a biological gene and is excluded from all biological denominators and subtype sets, while retained in the mapping audit. There are {int(eco_counts['vertices_gene.csv']):,} exported gene vertices and **{n_eligible:,} eligible gene symbols** (including symbols unobserved in a particular subtype).

### Supplementary recovery, not silent export repair

`ECoGraph/data/coad_mut_data.csv` contains {recovery.get('samples', 'NA')} distinct patient records and {recovery.get('genes', 'NA')} binary gene columns; patient IDs equal the first 12 characters of tumor barcodes. Its patient IDs match all 404 sample vertices. After excluding `Tumor_Sample_Barcode` and `patient_id`, gene columns match eligible exported symbols exactly. Values were verified to be only 0/1. Every exported subtype pair agrees with `coad_subtype_data.csv`; two additional GI.CIN assignments in that source lack corresponding exported samples and are not introduced.

The supplementary analysis joins the binary matrix to the **exported** 339 subtype assignments, counts each patient–gene pair once, and interprets 1 as an observed binary mutation flag according to the local processing documentation. It does not count variants, infer activating/inactivating status, or regard a 0 as proof that the gene was comprehensively tested. {recovery.get('observations_all_samples', 'NA')} observed pairs occur across all 404 matrix patients. No reconstruction is written into the source export. Export-only subtype molecular metrics remain **NA**, not zero, in `neptune_only_subtype_status.csv`.

## PrimeKG schema, identifiers, and relevant relationships

The raw CSV scan finds {prime_info['raw_rows']:,} ordered rows, {len(nodes):,} node indices, {nodes.type.nunique()} node types and {prime_info['relation_labels']} relation labels. Node types are {', '.join(sorted(nodes.type.unique()))}. Each edge row has `relation`, `display_relation` and endpoint `index`, `id`, `type`, `name`, `source` fields. All identifier fields are read as strings. No blank cells or conflicting metadata per node index were found. There are {len(dup_nodes)} node rows sharing `(type, source, id)` ({int(dup_nodes.type.eq('gene/protein').sum())} gene rows, {int(dup_nodes.type.eq('anatomy').sum())} anatomy rows); exact-ID aliases are counted once for biological matching, without modifying original nodes.

Gene/protein nodes use NCBI IDs ({int(prime_counts['gene/protein']):,} original nodes, {int(overall.unique_primekg_genes.iloc[0]):,} distinct IDs). Pathways use Reactome stable `R-HSA-...` identifiers ({int(prime_counts['pathway']):,}). Drugs use DrugBank `DB...` IDs ({int(prime_counts['drug']):,}). Diseases use MONDO and MONDO_grouped IDs ({int(prime_counts['disease']):,} nodes). `primekg_node_types.csv` and `primekg_relations.csv` contain the complete observed namespace and typed relation inventory.

Crucially, `drug_protein` has four display labels: **target, enzyme, carrier, transporter**. Only `display_relation == target` is counted as direct drug–target connectivity. Enzymes, carriers, and transporters are excluded from target counts. Drug→Gene and Gene→Drug rows are reciprocal representations, verified and deduplicated. Likewise, `indication` rows are deduplicated by DrugBank drug and disease ID; `contraindication` and `off-label use` are kept distinct and are not therapeutic positives in these counts.

{table(pd.read_csv(out / 'drug_gene_relation_counts.csv'))}

There are **{context['direct_target_pairs']:,} direct target pairs**, involving **{context['target_genes']:,} genes and {context['target_drugs']:,} drugs**. Context relations also include `protein_protein`, `disease_protein`, `pathway_protein`, GO functional annotations, anatomy expression, and positive/negative disease–phenotype links. These offer possible message-passing context, not proof of drug response.

## Gene identifier compatibility

ECoGraph lacks stable gene identifiers. Mapping uses exact, case-sensitive gene symbols against PrimeKG gene/protein `name`, accepting a match only if it resolves to one NCBI ID. Only leading/trailing whitespace may be normalized, with a separate status. We do not uppercase everything, remove punctuation, use fuzzy names, or resolve aliases automatically. Multiple indices with the same NCBI ID are one biological candidate; a symbol mapping to multiple NCBI IDs is ambiguous and is not accepted. Accepted symbol matches are provisional identifier bridges, not an authoritative HGNC crosswalk. A curated stable-ID crosswalk is recommended before integration.

{table(overall)}

The main percentage excludes the known `patient_id` artifact; the CSV also reports raw-export coverage so the denominator is explicit. `identifier_mapping_audit.csv` and `gene_mapping_audit.csv` retain exact, normalized, ambiguous, unmatched, and excluded rows, accepted IDs, and all candidate PrimeKG indices. A high overlap here establishes useful bridges, not subtype-specific therapeutic relevance.

## Subtype molecular associations and mapping

**All tables in this section are supplementary matrix-derived.** “Associated” means observed in at least one assigned sample, not enriched, causal, a driver, or exclusive biology. Sample–gene observations are unique binary patient–gene pairs.

{table(sub[["subtype", "samples", "associated_genes", "sample_gene_observations", "mapped_genes", "mapping_pct", "unmatched_genes", "ambiguous_genes", "exclusive_observed_genes"]])}

Per-gene counts and frequencies are provided in the four `subtype_gene_sets/*.csv` files. Frequency is the number of subtype patients with an observed gene flag divided by subtype sample count. Set size and exclusivity are highly sensitive to both sample count and mutation burden.

Selected concrete gene examples (descriptive frequencies, not drug-response biomarkers):

{table(example_frame)}

For example, APC can map into PrimeKG even though no direct drug-target edge is recorded for it; a recorded target connection for TP53 or KRAS does not establish that any drug is appropriate for the observed variant. Missing target edges are also not proof that a gene is untargetable.

### Mutation burden and enrichment screens

{table(burden)}

These burdens count distinct flagged genes per sample, **not mutations per megabase**. Hypermutated groups have more opportunities to accumulate passenger flags; GI.HM-SNV has only six assigned samples. We therefore do not treat a union of observed genes as a subtype-specific molecular signature.

We executed two explicitly exploratory one-sided tests for increased binary prevalence, comparing each subtype with all other typed patients: exact Fisher/hypergeometric tests, and a burden-stratified Cochran–Mantel–Haenszel normal-approximation screen over pooled burden quintiles. The latter sums observed-minus-expected counts and conditional variances in strata containing both groups; strata without overlap contribute no information. We suppress the burden screen unless each group has at least 20 samples and at least half of each group's samples fall in mixed strata. These are minimal operational guardrails, not proof of test validity. Disabled tests receive p=1 before BH; their summary results are NA. Unfiltered approximate p-values are retained only for diagnostics. BH correction is applied across all 4 × 21,332 hypotheses separately for each test family. The joint screen requires both q-values <0.05 and greater observed subtype prevalence. All tests and support counts are exported; this is not a validated enrichment procedure.

{table(enrichment)}

The HM-indel burden strata contain only eight comparison patients, so that burden screen is disabled; GI.HM-SNV also fails the sample-size and comparison-support guardrails. Naive Fisher testing produces thousands of hypermutated-group hits, illustrating why unadjusted prevalence must not be presented as subtype-specific biology. The CMH approximation can still be poor for sparse genes; burden quintiles provide incomplete adjustment and include each tested gene in the total burden. No callable-region, gene-length, variant-consequence, sequencing-depth, treatment, or external replication data are used. Thus any eligible joint candidates are **exploratory associations, not validated subtype-enriched genes**. Required validated-enrichment counts remain NA. Their mapped/target counts are separately labeled exploratory, rather than promoted to reliable enriched-gene coverage.

## Pathway compatibility

ECoGraph's {len(path_audit):,} endpoint IDs are collection-prefixed names: {collection_counts.get('REACTOME', 0):,} REACTOME, {collection_counts.get('BIOCARTA', 0):,} BIOCARTA, and {collection_counts.get('KEGG', 0):,} KEGG. Example: `pathway:BIOCARTA_41BB_PATHWAY`. PrimeKG has {int(prime_counts['pathway']):,} Reactome stable IDs such as `R-HSA-109581` (Apoptosis). ECoGraph does not supply corresponding stable database IDs or a name-to-ID crosswalk.

{table(pathway_summary.drop(columns=['interpretation']))}

No pathway is accepted solely on a similar or normalized name. The audit suggests candidate Reactome label correspondences by stripping `REACTOME_` and normalizing punctuation/case; these are explicitly unverified and never used to merge entities or claim overlap. **Zero confidently established mappings is not evidence of zero biological overlap. Actual pathway overlap percentage is not reliably computable from supplied identifiers.** BIOCARTA and KEGG pathways must remain separate from Reactome entities unless a defensible equivalence relation is supplied. A future crosswalk should preserve source collection, version, stable ID, and mapping provenance.

## Direct drug-target connectivity by subtype

{table(drug[["subtype", "mapped_genes", "mapped_genes_with_direct_targets", "direct_target_pct_of_mapped_genes", "unique_connected_drugs", "unique_direct_drug_gene_pairs"]])}

These are direct `drug_protein`/`target` edges, not indirect graph paths. Drug counts are distinct DrugBank IDs; gene coverage counts distinct ECoGraph symbols resolving to unambiguous PrimeKG IDs. The CSV separately records distinct mapped NCBI IDs. A connected drug may target an observed gene without modulating the relevant variant, acting in the needed direction, reaching the tissue, or improving disease. Connectivity, indication, prediction, and actual efficacy are separate quantities. There are no observed drug–GI subtype therapeutic labels in these inputs.

## Cross-subtype comparison

{table(overlap[["subtype_a", "subtype_b", "genes_a", "genes_b", "overlap", "jaccard"]])}

**{len(allshared):,} observed genes are shared across all four subtypes.** Exclusives are listed by count above; gene-set CSVs allow identification. `subtype_exploratory_screen_overlap.csv` repeats comparisons for joint screen candidates, without treating them as validated enrichment. Different observed unions and frequencies can give different input representations, but may primarily encode mutation burden and sample size. They do not establish different therapeutic vulnerabilities. Weighted prevalence representations and burden-controlled validation should precede any claim of meaningful subtype drug-ranking differences. Mapping percentages and direct-target percentages also differ, as quantified above.

## Therapeutic supervision and colorectal disease context

PrimeKG contains **{context['indication_pairs']:,} unique known indication pairs**, spanning **{context['indication_drugs']:,} drugs and {context['indication_diseases']:,} diseases**. These provide a technically usable supervised drug–disease relation, although abundance alone does not establish adequacy for the specific CRC subtypes. Unknown pairs are not verified negatives. A future evaluation must remove target/reverse/alias equivalents and prevent indication text or related relations from trivially revealing the answer; disease-family or source/temporal holdouts are preferable to a naive row split.

A lexical scan for `colorectal`, whole-word `colon`, or `rectal` yields {context['crc_lexical_candidates']} disease candidates, including colorectal cancer, colon adenocarcinoma, colorectal adenocarcinoma, rectal cancer, benign lesions, hereditary categories, and possible irrelevant name matches. This is an inclusive audit list, **not a curated equivalence mapping**. Across these candidates there are {context['crc_indication_pairs']} indication pairs and {context['crc_indication_drugs']} drugs; see `primekg_colorectal_disease_candidates.csv` and `primekg_crc_indications.csv` for exact records. GI.CIN, GI.HM-indel, GI.HM-SNV and GI.GS must not be merged into those broad disease nodes.

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

Run `python3 analysis/run_compatibility.py` from the project root. `analysis/README.md` documents dependencies and optional suppression of supplementary recovery. The script never modifies datasets or existing preprocessing outputs. PrimeKG is processed in chunks with only relevant drug edges retained; no second full edge list or merged graph is saved. All original mandatory inputs are SHA-256 checked before and after execution. `analysis/results/README.md` indexes artifacts. Statistical screening results are supplementary, and no drugs are recommended as treatments.
"""
    (out / "compatibility_report.md").write_text(report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primekg", type=Path, default=Path("kg.csv"))
    parser.add_argument("--ecograph", type=Path, default=Path("ECoGraph/data/neptune_csv"))
    parser.add_argument("--output", type=Path, default=Path("analysis/results"))
    parser.add_argument("--no-supplementary-recovery", action="store_true")
    args = parser.parse_args()
    out = args.output
    source_paths = [args.primekg] + [args.ecograph / f for f in FILES]
    if any(out.resolve() == p.resolve() or out.resolve() in p.resolve().parents for p in source_paths):
        raise ValueError("Outputs must not contain source files")
    out.mkdir(parents=True, exist_ok=True)
    manifest = inspect_sources(source_paths, out)
    print((out / "source_inspection.txt").read_text(), flush=True)
    frames = {name: read_csv(args.ecograph / name) for name in FILES}
    if len(frames["edges_has_gene_part001.csv"]):
        raise ValueError("The sample-gene export is now populated. Reinspect its value/schema semantics before using the blank-export recovery workflow.")
    eco_integrity = ecograph_audit(frames, out)
    nodes, edges, inventory = scan_primekg(args.primekg, out)
    prime_genes = nodes[nodes.type.eq("gene/protein")]
    assert set(prime_genes.source) == {"NCBI"}
    mapping = map_genes(frames["vertices_gene.csv"], prime_genes)
    mapping.to_csv(out / "gene_mapping_audit.csv", index=False)
    eligible = mapping.status.ne("excluded_metadata_artifact")
    matched = mapping.primekg_id.ne("")
    overall = pd.DataFrame([{"raw_export_gene_vertices": len(mapping), "eligible_ecograph_genes": int(eligible.sum()),
        "primekg_gene_nodes": len(prime_genes), "unique_primekg_genes": prime_genes.id.nunique(),
        "matching_ecograph_genes": int(matched.sum()), "matching_unique_ncbi_ids": mapping.loc[matched, "primekg_id"].nunique(),
        "ecograph_mapping_pct": pct(matched.sum(), eligible.sum()),
        "raw_export_mapping_pct": pct(matched.sum(), len(mapping)),
        "primekg_in_ecograph_pct": pct(mapping.loc[matched, "primekg_id"].nunique(), prime_genes.id.nunique()),
        "ambiguous_ecograph_genes": int(mapping.status.eq("ambiguous").sum()),
        "unmatched_eligible_genes": int(mapping.status.eq("unmatched").sum())}])
    overall.to_csv(out / "gene_overlap_summary.csv", index=False)
    pathway_mapping, pathway_summary = pathway_audit(frames["edges_in_pathway.csv"], nodes[nodes.type.eq("pathway")], out)
    pd.concat([mapping, pathway_mapping], ignore_index=True).to_csv(out / "identifier_mapping_audit.csv", index=False)
    targets, context, crc = therapeutic_context(edges, nodes, out)
    (out / "therapeutic_context_summary.json").write_text(json.dumps(context, indent=2, default=int) + "\n")
    x, labels, recovery = (None, None, {"status": "disabled"}) if args.no_supplementary_recovery else supplementary_matrix(args.ecograph, frames, mapping, out)
    (out / "supplementary_recovery_audit.json").write_text(json.dumps(recovery, indent=2) + "\n")
    sets, sub, drug, overlap, burden = subtype_analysis(frames, mapping, targets, x, labels, out)
    if x is not None:
        screens, enrichment = exploratory_enrichment(x, labels, mapping, targets, out)
        screened_sets = {s: set(screens.loc[screens.subtype.eq(s) & screens.screen_candidate, "symbol"]) for s in SUBTYPES}
        screen_overlap = compare_sets(screened_sets, "exploratory_screen_not_validated_enrichment")
        disabled = set(enrichment.loc[~enrichment.burden_screen_eligible, "subtype"])
        unavailable = screen_overlap.subtype_a.isin(disabled) | screen_overlap.subtype_b.isin(disabled)
        screen_overlap.loc[unavailable, ["genes_a", "genes_b", "overlap", "union", "jaccard"]] = np.nan
        screen_overlap["status"] = np.where(unavailable, "unavailable: burden-screen guardrails failed", "exploratory only")
        screen_overlap.to_csv(out / "subtype_exploratory_screen_overlap.csv", index=False)
    else:
        enrichment = pd.DataFrame([{"status": "not computable without sample-gene observations"}])
    build_report(out, eco_integrity, nodes, mapping, overall, pathway_summary,
                 sub, drug, overlap, burden, enrichment, context, recovery, crc)
    for item in manifest:
        assert sha256(Path(item["path"])) == item["sha256"], f"Source changed: {item['path']}"
    if recovery.get("status") == "verified_supplementary_reconstruction":
        assert sha256(args.ecograph.parent / "coad_mut_data.csv") == recovery["matrix_sha256"]
        assert sha256(args.ecograph.parent / "coad_subtype_data.csv") == recovery["subtype_sha256"]
    (out / "validation.json").write_text(json.dumps({"mandatory_source_hashes_unchanged": True,
        "supplementary_source_hashes_unchanged": True if x is not None else None,
        "pandas": pd.__version__, "numpy": np.__version__, "scipy": scipy.__version__, "verdict": "VIABLE WITH MODIFICATIONS",
        "export_only_subtype_molecular_status": "unavailable: sample-gene edge file empty",
        "supplementary_recovery": recovery["status"]}, indent=2) + "\n")
    files = sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file())
    (out / "README.md").write_text("# Compatibility analysis outputs\n\nStart with [compatibility_report.md](compatibility_report.md).\n\n"
        "Subtype metrics and gene sets are supplementary matrix-derived; Neptune-only results are unavailable.\n"
        "Blank numeric values mean unavailable/undefined, not zero. All enrichment screens are exploratory.\n\n"
        + "\n".join(f"- `{p}`" for p in files) + "\n")
    print("\nOVERALL GENE COVERAGE\n", overall.to_string(index=False))
    print("\nSUPPLEMENTARY SUBTYPE GENE COVERAGE\n", sub.to_string(index=False))
    print("\nSUPPLEMENTARY DIRECT DRUG-TARGET COVERAGE\n", drug.to_string(index=False))
    print("\nPATHWAYS\n", pathway_summary.to_string(index=False))
    print("\nCROSS-SUBTYPE OVERLAP\n", overlap.to_string(index=False))
    print("\nVERDICT: VIABLE WITH MODIFICATIONS; export repair required; prediction unvalidated")
    print("\nOUTPUTS\n" + "\n".join(str(p) for p in sorted(out.rglob("*")) if p.is_file()))
