"""Profile raw rows and export nodes and a phenotype/disease/gene subgraph.

Identifiers are opaque strings. No biological equivalence or directionality is
inferred, and raw endpoint order is preserved.
"""

import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd

NODE_FIELDS = ["index", "id", "type", "name", "source"]
EDGE_FIELDS = ["x_index", "y_index", "relation", "display_relation"]
COLUMNS = ["relation", "display_relation"] + [
    f"{side}_{field}" for side in ("x", "y") for field in NODE_FIELDS
]
CORE_TYPES = {"effect/phenotype", "disease", "gene/protein"}
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def counter_frame(counter, columns, count_name="raw_row_count"):
    return pd.DataFrame(
        [(*key, value) for key, value in sorted(counter.items())],
        columns=[*columns, count_name],
    )


def prepare(input_path, output_dir, report_dir, chunksize=250_000):
    input_path, output_dir, report_dir = map(Path, (input_path, output_dir, report_dir))
    if input_path.resolve().parent == output_dir.resolve():
        raise ValueError("Use a separate output directory for derived data.")
    header = pd.read_csv(input_path, nrows=0).columns.tolist()
    if header != COLUMNS:
        raise ValueError(f"Unexpected CSV schema: {header}; expected {COLUMNS}")
    output_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    node_parts, core_parts = [], []
    relationships, missing = Counter(), Counter()
    rows = self_loops = core_raw_rows = 0
    # Disable default NA tokens: names/IDs such as 'NA' must stay literal.
    for chunk in pd.read_csv(input_path, dtype="string", keep_default_na=False,
                             chunksize=chunksize):
        rows += len(chunk)
        missing.update(chunk.eq("").sum().to_dict())
        for field in ("x_index", "y_index"):
            chunk[field] = pd.to_numeric(chunk[field], errors="raise").astype("int64")
        self_loops += int(chunk.x_index.eq(chunk.y_index).sum())
        relationships.update(chunk.groupby(
            ["relation", "display_relation", "x_type", "y_type"], dropna=False
        ).size().to_dict())
        for side in ("x", "y"):
            nodes = chunk[[f"{side}_{f}" for f in NODE_FIELDS]].copy()
            nodes.columns = NODE_FIELDS
            node_parts.append(nodes.drop_duplicates())
        mask = chunk.x_type.isin(CORE_TYPES) & chunk.y_type.isin(CORE_TYPES)
        core = chunk.loc[mask, EDGE_FIELDS]
        core_raw_rows += len(core)
        core_parts.append(core.groupby(EDGE_FIELDS, dropna=False).size()
                          .rename("raw_row_count").reset_index())

    observations = pd.concat(node_parts, ignore_index=True).drop_duplicates()
    observations = observations.sort_values(NODE_FIELDS).reset_index(drop=True)
    conflicts = observations[observations.duplicated("index", keep=False)]
    conflicts.to_csv(report_dir / "node_metadata_conflicts.csv", index=False)
    if len(conflicts):
        raise ValueError("Conflicting metadata for a node index; inspect node_metadata_conflicts.csv")
    nodes = observations
    nodes.to_csv(output_dir / "nodes.csv", index=False)
    core_nodes = nodes[nodes.type.isin(CORE_TYPES)]
    core_nodes.to_csv(output_dir / "phenotype_disease_gene_nodes.csv", index=False)
    edges = pd.concat(core_parts, ignore_index=True).groupby(
        EDGE_FIELDS, as_index=False, dropna=False
    ).raw_row_count.sum().sort_values(EDGE_FIELDS).reset_index(drop=True)
    edges.to_csv(output_dir / "phenotype_disease_gene_edges.csv", index=False)
    node_ids = set(nodes["index"])
    assert set(edges.x_index).union(edges.y_index) <= node_ids
    assert int(edges.raw_row_count.sum()) == core_raw_rows
    assert len(nodes) == nodes["index"].nunique()
    relation_table = counter_frame(relationships,
        ["relation", "display_relation", "x_type", "y_type"])
    relation_table.to_csv(report_dir / "relationships.csv", index=False)
    node_counts = nodes.groupby("type").size().rename("node_count").reset_index()
    node_counts.to_csv(report_dir / "node_types.csv", index=False)
    aliases = nodes[nodes.duplicated(["type", "source", "id"], keep=False)]
    aliases.to_csv(report_dir / "repeated_external_ids.csv", index=False)

    # Reciprocal rows are reported, not collapsed: semantics vary by relation.
    reverse = edges.rename(columns={"x_index": "y_index", "y_index": "x_index"})
    paired = edges.merge(reverse[EDGE_FIELDS], on=EDGE_FIELDS, how="inner")
    summary = {
        "input": str(input_path), "input_bytes": input_path.stat().st_size,
        "columns": header, "raw_rows": rows, "unique_nodes": len(nodes),
        "node_types": len(node_counts), "relation_labels": len({k[0] for k in relationships}),
        "typed_relationships": len(relation_table), "blank_cells_by_column": dict(missing),
        "raw_self_loop_rows": self_loops, "node_metadata_conflict_rows": len(conflicts),
        "repeated_external_id_node_rows": len(aliases),
        "core_types": sorted(CORE_TYPES), "core_nodes": len(core_nodes),
        "core_raw_rows": core_raw_rows, "core_unique_ordered_edges": len(edges),
        "core_duplicate_ordered_rows": core_raw_rows - len(edges),
        "core_edges_with_reverse_same_labels": len(paired),
        "full_graph_duplicate_rows": "not computed; raw relationship counts include duplicates",
    }
    (report_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = ["# Supplied PrimeKG dataset profile", "",
        f"Raw rows: {rows:,}. Unique nodes: {len(nodes):,}.", "",
        "## Node types", "", "| Type | Nodes |", "| --- | ---: |"]
    lines.extend(f"| {r.type} | {r.node_count:,} |" for r in node_counts.itertuples())
    lines.extend(["", "## Typed relationships (raw row counts)", "",
        "| Relation | Display label | X type | Y type | Rows |", "| --- | --- | --- | --- | ---: |"])
    lines.extend(f"| {r.relation} | {r.display_relation} | {r.x_type} | {r.y_type} | {r.raw_row_count:,} |"
                 for r in relation_table.itertuples())
    lines.extend(["", "## Phenotype–disease–gene extraction", "",
        f"{len(core_nodes):,} nodes; {core_raw_rows:,} raw rows; {len(edges):,} unique ordered edges.",
        f"Removed {core_raw_rows-len(edges):,} repeated ordered rows; multiplicity is retained in `raw_row_count`.",
        f"{len(paired):,} exported ordered edges have a reverse edge with identical relation labels.",
        "", "See ../../docs/research/primekg_exploration.md for semantics, limitations, and construction examples.", ""])
    (report_dir / "profile.md").write_text("\n".join(lines))
    print(json.dumps(summary, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=PROJECT_ROOT / "data/raw/primekg/kg.csv")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "data/processed/primekg")
    parser.add_argument("--report-dir", type=Path, default=PROJECT_ROOT / "reports/primekg")
    parser.add_argument("--chunksize", type=int, default=250_000)
    args = parser.parse_args()
    if args.chunksize <= 0:
        parser.error("--chunksize must be positive")
    prepare(args.input, args.output_dir, args.report_dir, args.chunksize)
