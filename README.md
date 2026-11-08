# PrimeKG exploration

A small pandas workflow for understanding the supplied `kg.csv` and preparing
phenotype–disease–gene data. It does not implement diagnosis, a GNN, or an agent.
All findings refer to this local CSV; no external ontology or clinical evidence
has been added.

For the follow-up coverage analysis, masked-edge baselines, limitations, and
project recommendation, read [primekg_findings.md](primekg_findings.md).
Reproduce those results after preprocessing with
`python3 scripts/explore_inference.py`; small outputs go to `reports/inference/`.

## Findings from the supplied file

- 8,100,498 raw edge rows, 129,375 unique node indices, 10 node types,
  and 30 relation labels. There are 56 combinations of relation, display label,
  and endpoint types. In particular, `drug_protein` has multiple display labels
  (`target`, `enzyme`, `carrier`, `transporter`), so display labels carry information.
- The core contains 17,080 disease nodes (15,813 `MONDO`, 1,267
  `MONDO_grouped`), 15,311 phenotype nodes (`HPO`), and 27,671 gene/protein
  nodes (`NCBI`). Grouped disease IDs are preserved as opaque strings.
- Disease-to-phenotype rows include 150,317 positive and 1,193 negative
  links. Disease-to-gene rows number 80,411; phenotype-to-gene rows number
  3,330. Each also has a reciprocal row. These are dataset associations,
  not necessarily distinct diseases or independently supported evidence items.
- All 1,214,512 selected ordered edges are unique and have a reciprocal
  edge with identical labels. This corresponds to 607,256 unordered labeled
  pairs, without establishing that every relation is biologically undirected.
- No blank cells, raw self-loops, or conflicting metadata per node index were
  found. There are 126 node rows sharing `(type, source, id)` with another
  index; inspect `reports/repeated_external_ids.csv` before any ID-based merging.

The full node and relation inventory is in [reports/profile.md](reports/profile.md).

## Run

From the repository root, with Python 3.10 or later:

```sh
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 scripts/prepare_primekg.py
```

The script reads in chunks (default 250,000 rows), then consolidates node metadata
and the selected edges in memory. It is not a constant-memory pipeline. Use
`--chunksize`, `--input`, `--output-dir`, and `--report-dir` to customize it.
The raw CSV is read only. Raw and derived large CSVs are excluded from Git.

## CSV structure and graph identity

Each row represents a labeled edge from the `x` endpoint to the `y` endpoint.
There are 12 columns:

| Columns | Meaning and handling |
| --- | --- |
| `relation`, `display_relation` | Machine relation label and human-facing label; both retained |
| `x_index`, `y_index` | Integer node keys used to join endpoints to the node table |
| `x_id`, `y_id` | External identifiers; read as strings without numeric normalization |
| `x_type`, `y_type` | Endpoint node types |
| `x_name`, `y_name` | Display names, not unique identifiers |
| `x_source`, `y_source` | Node identifier namespaces/sources; these are not per-edge evidence citations |

`nodes.csv` consolidates both endpoint columns into one row per `index` with
`id`, `type`, `name`, and `source`. The script checks metadata consistency across
every occurrence of each node. Conflicts are exported and stop processing;
they are never silently resolved. It also reports repeated `(type, source, id)`
values across distinct indices instead of merging nodes by name or external ID.

Use the supplied index as graph identity within this dataset. Do not assume
indices are contiguous or stable across other dataset versions. If a future
library needs consecutive indices, create and retain an explicit mapping.
Names and IDs are preserved exactly, including leading zeros and literal `NA`.
Blank fields are counted; metadata is not filled or dropped.

## Outputs

| File | Contents |
| --- | --- |
| `data/processed/nodes.csv` | All graph nodes and original metadata |
| `data/processed/phenotype_disease_gene_nodes.csv` | All nodes of the three selected types, including nodes with no selected edges |
| `data/processed/phenotype_disease_gene_edges.csv` | Edges whose endpoints both have selected types; unique ordered `(x_index, y_index, relation, display_relation)` tuples with `raw_row_count` |
| `reports/profile.md` | Dataset counts, node types, and typed relationships |
| `reports/summary.json` | Machine-readable counts, blanks, and validation findings |
| `reports/relationships.csv` | Raw row counts by relation, display label, and endpoint types |
| `reports/node_types.csv` | Unique node counts by type |
| `reports/node_metadata_conflicts.csv` | Conflicting metadata observations, empty when consistent |
| `reports/repeated_external_ids.csv` | Distinct node indices sharing a typed, namespaced external ID |

The complete graph's edges remain available in `kg.csv`; they are not copied to
another large file. Full-graph relationship counts include all raw rows and do
not claim to measure unique biological associations. Exact ordered duplicates
are consolidated only in the selected subgraph, retaining their multiplicity.
`raw_row_count` is a storage count, not a biological confidence or evidence score.

## Relevant relationships and interpretation

The extraction includes **all** relations between `effect/phenotype`, `disease`,
and `gene/protein`, leaving later scope decisions open:

- `disease_phenotype_positive`: phenotype associations for exploring candidate diseases.
- `disease_phenotype_negative`: explicitly negative associations, kept separate.
- `disease_protein`: disease–gene/protein associations for candidate gene exploration.
- `phenotype_protein`: direct phenotype–gene/protein associations.
- `protein_protein`, `phenotype_phenotype`, `disease_disease`: possible surrounding graph context.

These labels alone do not establish causality, inheritance, diagnostic
sensitivity, prevalence, or the precise semantics of within-type edges.
Positive and negative dataset relations are knowledge-graph annotations, not
patient observations. A symptom that was not entered by a patient is unknown,
not an explicitly absent finding. Missing edges are also not negative evidence.

The CSV includes reciprocal rows for many associations. Preserve their supplied
orientation initially; do not add reverse edges blindly or treat reverse copies
as independent evidence. Same-type relations may have direction-specific
semantics requiring source documentation before changing their direction.
If a later experiment splits edges for evaluation, keep reciprocal copies of
the same association in one split to avoid leakage.

Drug, pathway, anatomy, exposure, and functional ontology data remain represented
in the full node table and relationship inventory. They may be useful later,
but are not required for a first phenotype–disease–gene exploration. The CSV
does not supply patient variants, patient cases, a dedicated rare-disease flag,
or a clinical test catalog. Candidate ranking and test recommendations would
need additional design and evidence beyond this preprocessing.

## Explore with pandas

```python
import pandas as pd

nodes = pd.read_csv("data/processed/nodes.csv", dtype={
    "index": "int64", "id": "string", "type": "string",
    "name": "string", "source": "string",
}, keep_default_na=False)
edges = pd.read_csv("data/processed/phenotype_disease_gene_edges.csv",
                    keep_default_na=False)

# Resolve a phenotype interactively: a name search is not a clinical mapping.
matches = nodes[nodes["type"].eq("effect/phenotype") &
                nodes["name"].str.contains("seizure", case=False, regex=False)]
print(matches[["index", "id", "name", "source"]].head(20))

# Inspect all positive disease associations touching the matched nodes.
links = edges[edges.relation.eq("disease_phenotype_positive") &
              (edges.x_index.isin(matches["index"]) |
               edges.y_index.isin(matches["index"]))]
named = links.merge(nodes.add_prefix("x_"), on="x_index", validate="many_to_one")
named = named.merge(nodes.add_prefix("y_"), on="y_index", validate="many_to_one")
print(named[["x_name", "x_type", "relation", "y_name", "y_type"]].head())
```

The node table plus labeled endpoint tuples are sufficient to construct a
heterogeneous graph later: group edges by `(x_type, relation, y_type)` after
joining endpoint metadata. Keep multiple relation labels between the same pair;
a plain unlabeled adjacency matrix would discard information. No graph library
is needed for the current exploration.
