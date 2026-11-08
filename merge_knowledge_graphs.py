#!/usr/bin/env python3
"""Merge two CSV knowledge graphs with conservative identifier matching.

Supported base graph:
  --base-format primekg: PrimeKG kg.csv with x_/y_ endpoint columns.
  --base-format tables: nodes.csv and edges.csv (see README).

Incoming graph: nodes.csv and relationships.csv (Neo4j/SBML-style columns).
Requires pandas. No network access, ontology guessing, or inferred edges.
"""
from __future__ import annotations
import argparse
import ast
import json
import re
from pathlib import Path
import pandas as pd


def read(path, required=True):
    if not path or not Path(path).is_file():
        if required:
            raise FileNotFoundError(f"Missing CSV: {path}")
        return pd.DataFrame()
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def column(df, names, default=""):
    for name in names:
        if name in df.columns:
            return df[name].fillna("").astype(str)
    return pd.Series([default] * len(df), index=df.index, dtype=str)


def parse(value):
    if isinstance(value, dict):
        return value
    if not value:
        return {}
    for f in (json.loads, ast.literal_eval):
        try:
            result = f(value)
            return result if isinstance(result, dict) else {}
        except (ValueError, SyntaxError, TypeError):
            continue
    return {}


def standard_nodes(df, graph):
    ids = column(df, ['node_id', 'id', 'element_id', '~id'])
    if not ids.any():
        raise ValueError(f'{graph} nodes require node_id or id')
    props = column(df, ['properties']).map(parse)
    def value(row, keys):
        for k in keys:
            if k in df.columns and str(row[k]).strip():
                return str(row[k]).strip()
        lowered = {str(k).lower(): v for k, v in parse(row.get('properties', '')).items()}
        for k in keys:
            v = lowered.get(k.lower(), '')
            if v is not None and str(v).strip():
                return str(v).strip()
        return ''
    rows = []
    for _, r in df.iterrows():
        rows.append(dict(local_id=value(r, ['node_id', 'id', 'element_id', '~id']),
                         name=value(r, ['name', 'display_name', 'label', 'symbol']),
                         type=value(r, ['node_type', 'type', 'category', 'labels']),
                         identifier=value(r, ['identifier', 'curie', 'gene_id', 'entrez_id', 'ncbi_gene_id']),
                         source=value(r, ['source', 'database', 'namespace', 'identifier_source']),
                         symbol=value(r, ['symbol', 'gene_symbol']),
                         properties=json.dumps(parse(r.get('properties', '')), ensure_ascii=False)))
    out = pd.DataFrame(rows)
    if out.local_id.eq('').any() or out.local_id.duplicated().any():
        raise ValueError(f'{graph} nodes have empty or duplicate IDs')
    out['graph'] = graph
    return out


def standard_edges(df, graph):
    if df.empty and len(df.columns) == 0:
        raise ValueError(f'{graph} relationships CSV is missing headers')
    aliases = {'source_id': ['source_id', 'from_id', '~from', 'start_id', 'x_index'],
               'target_id': ['target_id', 'to_id', '~to', 'end_id', 'y_index'],
               'relation': ['relation', 'relationship', 'type', '~label']}
    for key, options in aliases.items():
        if not any(k in df.columns for k in options):
            raise ValueError(f'{graph} edges missing {key}; expected one of {options}')
    out = pd.DataFrame({k: column(df, opts) for k, opts in aliases.items()})
    out['properties'] = column(df, ['properties'])
    out['graph'] = graph
    out = out.loc[out.source_id.ne('') & out.target_id.ne('') & out.relation.ne('')].copy()
    return out


def primekg_graph(path):
    raw = read(path)
    required = ['x_index', 'y_index', 'relation']
    missing = [k for k in required if k not in raw]
    if missing:
        raise ValueError(f'PrimeKG file missing columns: {missing}')
    parts = []
    for side in ['x', 'y']:
        part = pd.DataFrame({'local_id': column(raw, [f'{side}_index']),
                             'name': column(raw, [f'{side}_name']),
                             'type': column(raw, [f'{side}_type']),
                             'identifier': column(raw, [f'{side}_id']),
                             'source': column(raw, [f'{side}_source']),
                             'symbol': '', 'properties': '{}', 'graph': 'base'})
        parts.append(part)
    nodes = pd.concat(parts, ignore_index=True)
    if nodes.local_id.eq('').any():
        raise ValueError('PrimeKG has empty endpoint indices')
    # One PrimeKG index must have consistent identity; fill blanks from other occurrences.
    for col in ['name', 'type', 'identifier', 'source']:
        conflicts = nodes.loc[nodes[col].ne('')].groupby('local_id')[col].nunique()
        if (conflicts > 1).any():
            raise ValueError(f'Conflicting PrimeKG metadata for {col} at indices {conflicts[conflicts > 1].index[:5].tolist()}')
    nodes = nodes.groupby('local_id', sort=False, as_index=False).agg(
        {c: (lambda s: next((x for x in s if x), '')) for c in ['name','type','identifier','source','symbol']}
        | {'properties': 'first', 'graph': 'first'})
    edges = standard_edges(raw, 'base')
    return nodes, edges


def identifier_key(source, identifier):
    s, i = str(source).strip().lower(), str(identifier).strip()
    if not i:
        return ''
    # Normalize widely used NCBI prefixes, but never equate unrelated databases.
    if i.lower().startswith(('ncbi:', 'ncbigene:', 'entrez:', 'geneid:')):
        i = i.split(':', 1)[1]
        s = 'ncbi'
    if s in ('ncbi', 'ncbigene', 'entrez', 'ncbi gene', 'ncbi_gene'):
        s = 'ncbi'
    if s == 'hgnc' and i.upper().startswith('HGNC:'):
        i = i.split(':', 1)[1]
    return f'{s}::{i.lower()}' if s else ''


def merge(base_nodes, base_edges, incoming_nodes, incoming_edges, explicit_map=None):
    base_nodes, incoming_nodes = base_nodes.copy(), incoming_nodes.copy()
    base_nodes['key'] = [identifier_key(s, i) for s, i in zip(base_nodes.source, base_nodes.identifier)]
    incoming_nodes['key'] = [identifier_key(s, i) for s, i in zip(incoming_nodes.source, incoming_nodes.identifier)]
    # Only unique base keys are eligible. Never match solely by display name.
    unique = base_nodes.loc[base_nodes.key.ne('') & ~base_nodes.key.duplicated(keep=False)]
    by_key = dict(zip(unique.key, unique.local_id))
    by_id = set(base_nodes.local_id)
    overrides = {}
    if explicit_map is not None:
        for r in explicit_map.itertuples(index=False):
            if r.incoming_id in overrides:
                raise ValueError(f'Duplicate mapping override: {r.incoming_id}')
            if r.base_id not in by_id:
                raise ValueError(f'Mapping override refers to unknown base_id: {r.base_id}')
            overrides[r.incoming_id] = r.base_id
    unknown = set(overrides) - set(incoming_nodes.local_id)
    if unknown:
        raise ValueError(f'Mapping overrides refer to unknown incoming IDs: {list(unknown)[:5]}')

    # Output identifiers are strings with prefixes; avoids collisions between input graphs.
    base_nodes['merged_id'] = 'base::' + base_nodes.local_id
    base_merged = dict(zip(base_nodes.local_id, base_nodes.merged_id))
    records = []
    for r in incoming_nodes.itertuples(index=False):
        if r.local_id in overrides:
            target, method = overrides[r.local_id], 'explicit_override'
        elif r.key and r.key in by_key:
            target, method = by_key[r.key], 'unique_source_identifier'
        else:
            target, method = '', 'new_node'
        records.append({'incoming_id': r.local_id, 'merged_id': base_merged[target] if target else 'incoming::' + r.local_id,
                        'base_id': target, 'mapping_method': method, 'identifier_key': r.key})
    audit = pd.DataFrame(records)
    incoming_lookup = dict(zip(audit.incoming_id, audit.merged_id))
    fresh = incoming_nodes.loc[incoming_nodes.local_id.isin(audit.loc[audit.mapping_method.eq('new_node'), 'incoming_id'])].copy()
    fresh['merged_id'] = 'incoming::' + fresh.local_id
    merged_nodes = pd.concat([base_nodes, fresh], ignore_index=True).drop(columns=['key'])
    merged_nodes = merged_nodes[['merged_id','graph','local_id','name','type','identifier','source','symbol','properties']]

    def remap(edges, lookup):
        e = edges.copy()
        e['merged_source_id'] = e.source_id.map(lookup)
        e['merged_target_id'] = e.target_id.map(lookup)
        return e
    base_remapped = remap(base_edges, base_merged)
    incoming_remapped = remap(incoming_edges, incoming_lookup)
    all_edges = pd.concat([base_remapped, incoming_remapped], ignore_index=True)
    rejected = all_edges.loc[all_edges.merged_source_id.isna() | all_edges.merged_target_id.isna()].copy()
    good = all_edges.dropna(subset=['merged_source_id','merged_target_id']).copy()
    # Keep graph provenance and properties; identical triples from different graphs remain traceable.
    good = good[['merged_source_id','relation','merged_target_id','graph','properties','source_id','target_id']]
    good = good.drop_duplicates(ignore_index=True)
    return merged_nodes, good, audit, rejected


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--base-format', choices=['primekg', 'tables'], default='primekg')
    p.add_argument('--base-edges', type=Path, required=True, help='PrimeKG kg.csv or generic base edges CSV')
    p.add_argument('--base-nodes', type=Path, help='Required for --base-format tables')
    p.add_argument('--incoming-nodes', type=Path, required=True)
    p.add_argument('--incoming-edges', type=Path, required=True)
    p.add_argument('--mapping', type=Path, help='Optional CSV: incoming_id,base_id (explicit curated matches)')
    p.add_argument('--out-dir', type=Path, default=Path('merged_graph'))
    args = p.parse_args()
    if args.base_format == 'primekg':
        bn, be = primekg_graph(args.base_edges)
    else:
        if not args.base_nodes:
            p.error('--base-nodes is required for --base-format tables')
        bn = standard_nodes(read(args.base_nodes), 'base')
        be = standard_edges(read(args.base_edges), 'base')
    an = standard_nodes(read(args.incoming_nodes), 'incoming')
    ae = standard_edges(read(args.incoming_edges), 'incoming')
    override = None
    if args.mapping:
        override = read(args.mapping)
        if not {'incoming_id', 'base_id'} <= set(override.columns):
            raise ValueError('Mapping CSV requires incoming_id,base_id')
    nodes, edges, audit, rejected = merge(bn, be, an, ae, override)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    nodes.to_csv(args.out_dir / 'nodes.csv', index=False)
    edges.to_csv(args.out_dir / 'edges.csv', index=False)
    audit.to_csv(args.out_dir / 'mapping_audit.csv', index=False)
    rejected.to_csv(args.out_dir / 'rejected_edges.csv', index=False)
    report = {'base_nodes': len(bn), 'base_edges': len(be), 'incoming_nodes': len(an),
              'incoming_edges': len(ae), 'matched_nodes': int(audit.mapping_method.ne('new_node').sum()),
              'new_nodes': int(audit.mapping_method.eq('new_node').sum()),
              'merged_nodes': len(nodes), 'merged_edges': len(edges), 'rejected_edges': len(rejected)}
    (args.out_dir / 'report.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print('Output:', args.out_dir.resolve())

if __name__ == '__main__':
    main()
