"""Read-only adapters for saved CRC experiments and derived PrimeKG tables.

Uses only the standard library: the CLI does not need torch, numpy or pandas.
All model results are cached outputs, never presented as fresh inference.
"""

import csv
import gzip
import json
import re
from functools import cached_property
from pathlib import Path

from .config import PROJECT_ROOT

RUN_ID = "20261108_034016"
SUBTYPES = ("GI.CIN", "GI.HM-indel", "GI.HM-SNV", "GI.GS")


class DataUnavailable(Exception):
    """A required local artifact is absent or not readable."""


def read_rows(path):
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))
    except (OSError, csv.Error, UnicodeError) as error:
        raise DataUnavailable(f"Cannot read {path.name}") from error


def normalize_subtype(value):
    value = value.strip().removeprefix("subtype:")
    return next((s for s in SUBTYPES if s.casefold() == value.casefold()), None)


class ResearchRepository:
    def __init__(self, root=PROJECT_ROOT):
        self.root = Path(root)
        # Pin a reviewed run rather than silently choosing a new, unreviewed model.
        self.run = self.root / "reports/crc_experiments" / RUN_ID

    def source(self, path):
        return str(path.relative_to(self.root))

    @cached_property
    def nodes(self):
        path = self.root / "data/processed/primekg/nodes.csv"
        rows = read_rows(path)
        if not rows or not {"index", "id", "name", "type", "source"} <= rows[0].keys():
            raise DataUnavailable("PrimeKG node table has an unexpected schema")
        return {row["index"]: row for row in rows}

    @cached_property
    def targets(self):
        return [r for r in read_rows(self.root / "reports/compatibility/primekg_direct_targets.csv")
                if r["display_relation"] == "target"]

    @cached_property
    def profiles(self):
        return read_rows(self.run / "profile_gene_evidence.csv.gz")

    def get_model_evaluation(self):
        path = self.run / "evaluation_metrics.csv"
        rows = read_rows(path)
        grouped = {}
        for row in rows:
            grouped.setdefault(row["condition"], []).append(float(row["test_ap"]))
        return {"status": "ok", "run_id": RUN_ID, "source": self.source(path),
            "task": "held-out PrimeKG indication recovery against sampled unknown pairs",
            "test_known_indications": 853, "test_unknown_pairs": 853,
            "metric": "average precision; not prediction accuracy or probability of efficacy",
            "conditions": {name: {"runs": len(values), "mean_test_ap": sum(values)/len(values),
                "min_test_ap": min(values), "max_test_ap": max(values)} for name, values in grouped.items()},
            "selected_cached_model": {"condition": "target_only", "seed": 42,
                "test_ap": next(float(r["test_ap"]) for r in rows if r["condition"] == "target_only" and r["seed"] == "42")},
            "interpretation": "Popularity outperformed the tested GNNs; molecular context gave no consistent AP gain over the drug-only ablation.",
            "subtype_efficacy_evaluated": False, "model_uses_phenotype_edges": False}

    def get_drug_predictions(self, subtype, top_k=10):
        subtype = normalize_subtype(subtype)
        if subtype is None:
            return {"status": "unsupported_subtype", "supported_subtypes": list(SUBTYPES)}
        if type(top_k) is not int or not 1 <= top_k <= 20:
            return {"status": "invalid_input", "message": "top_k must be an integer from 1 to 20"}
        path = self.run / "tie_aware_top_candidates.csv"
        selected = [r for r in read_rows(path) if r["seed"] == "42" and r["profile"] == "frequency"
                    and r["subtype"] == "subtype:" + subtype and int(r["tied_rank_min"]) <= top_k]
        if not selected:
            raise DataUnavailable("No matching reviewed cached prediction rows")
        selected.sort(key=lambda r: (int(r["tied_rank_min"]), r["drug_id"]))
        predictions = [{"drug_id": r["drug_id"], "drug_name": r["drug_name"],
            "rank_min": int(r["tied_rank_min"]), "tie_group_size": int(r["tie_group_size"]),
            "raw_logit": float(r["raw_logit"]),
            "recorded_target_ncbi_ids": r["recorded_target_ncbi_ids"].split("|"),
            "recorded_indication_for_crc_lexical_candidate": r["recorded_indication_for_crc_lexical_candidate"] == "True"}
            for r in selected]
        return {"status": "ok", "subtype": subtype, "backend": "cached_gnn_experiment",
            "fresh_inference": False, "run_id": RUN_ID, "model_condition": "target_only", "seed": 42,
            "profile": "frequency", "source": self.source(path), "requested_top_k": top_k,
            "returned_count": len(predictions), "ties_included": True, "predictions": predictions,
            "warning": "All candidates tied at the cutoff are returned, possibly more than top_k. Within-tie compound ordering is unresolved. Raw logits are not efficacy probabilities; graph nonmembership is not literature novelty.",
            "evaluation": self.get_model_evaluation()}

    def get_drug_targets(self, drug_id):
        drug_id = drug_id.strip()
        if not re.fullmatch(r"DB\d{5}", drug_id):
            return {"status": "invalid_input", "message": "Use a DrugBank ID such as DB00596"}
        drug = next((n for n in self.nodes.values() if n["type"] == "drug" and n["id"] == drug_id), None)
        if drug is None:
            return {"status": "unknown_drug", "drug_id": drug_id}
        gene_ids = sorted({r["gene_id"] for r in self.targets if r["drug_id"] == drug_id})
        genes = {}
        for row in self.nodes.values():
            if row["type"] == "gene/protein" and row["id"] in gene_ids:
                genes.setdefault(row["id"], {"gene_id": row["id"], "namespace": row["source"],
                    "symbols": set(), "node_indices": []})
                genes[row["id"]]["symbols"].add(row["name"])
                genes[row["id"]]["node_indices"].append(int(row["index"]))
        targets = [{**g, "symbols": sorted(g["symbols"])} for g in genes.values()]
        return {"status": "ok", "drug_id": drug_id, "drug_name": drug["name"],
            "relationship": {"relation": "drug_protein", "display_relation": "target"},
            "targets": targets, "target_count": len(targets),
            "source": "reports/compatibility/primekg_direct_targets.csv",
            "evidence_type": "recorded_direct_graph_relationship",
            "warning": "Excludes enzymes/carriers/transporters. Action direction, subtype benefit and clinical efficacy are not established."}

    def get_subtype_genes(self, subtype, limit=20, gene=""):
        subtype = normalize_subtype(subtype)
        if subtype is None:
            return {"status": "unsupported_subtype", "supported_subtypes": list(SUBTYPES)}
        if type(limit) is not int or not 1 <= limit <= 50:
            return {"status": "invalid_input", "message": "limit must be an integer from 1 to 50"}
        rows = [r for r in self.profiles if r["profile"] == "frequency" and r["subtype"] == "subtype:"+subtype]
        if not rows:
            raise DataUnavailable("No matching subtype profile rows")
        total_genes = len(rows)
        if gene.strip():
            gene = gene.strip()
            rows = [r for r in rows if (self.nodes[r["gene_node_index"]]["id"] == gene
                    if gene.isdigit() else r["symbol"] == gene)]
        rows.sort(key=lambda r: (-float(r["frequency"]), r["symbol"]))
        genes = []
        for row in rows[:limit]:
            node = self.nodes[row["gene_node_index"]]
            genes.append({"gene_id": node["id"], "namespace": node["source"], "symbol": row["symbol"],
                "gene_node_index": int(row["gene_node_index"]), "observed_samples": int(row["n_observed"]),
                "profiled_samples": int(row["n_profiled"]), "frequency": float(row["frequency"]),
                "pooling_weight": float(row["pooling_weight"])})
        return {"status": "ok", "subtype": subtype, "genes": genes, "total_genes": total_genes,
            "gene_filter": gene, "matching_genes": len(rows),
            "returned_count": len(genes), "truncated": len(rows) > limit, "profile": "frequency",
            "source": self.source(self.run / "profile_gene_evidence.csv.gz"),
            "evidence_type": "cohort_observation_not_validated_enrichment",
            "warning": "Observed flags are not driver, dependency or therapeutic-effect labels. Cohort sizes/burdens differ. No match means absent from the retained positive profile, not proof of biological absence."}

    @cached_property
    def gene_phenotypes(self):
        path = self.root / "data/processed/primekg/phenotype_disease_gene_edges.csv"
        result = {}
        try:
            with path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    if row["relation"] != "phenotype_protein":
                        continue
                    x, y = self.nodes[row["x_index"]], self.nodes[row["y_index"]]
                    gene, phenotype = (x, y) if x["type"] == "gene/protein" else (y, x)
                    if gene["type"] != "gene/protein" or phenotype["type"] != "effect/phenotype":
                        raise DataUnavailable("Unexpected phenotype-protein endpoint types")
                    result.setdefault(gene["id"], {})[phenotype["index"]] = phenotype
        except (OSError, KeyError, csv.Error) as error:
            raise DataUnavailable("Cannot read phenotype-protein graph edges") from error
        return result

    def get_gene_phenotypes(self, gene, limit=10):
        if type(limit) is not int or not 1 <= limit <= 50:
            return {"status": "invalid_input", "message": "limit must be an integer from 1 to 50"}
        gene = gene.strip()
        matches = [n for n in self.nodes.values() if n["type"] == "gene/protein" and
                   ((n["id"] == gene) if gene.isdigit() else (n["name"] == gene))]
        ids = sorted({n["id"] for n in matches})
        if not ids:
            return {"status": "unknown_gene", "query": gene, "message": "Use an NCBI ID or exact recorded symbol"}
        if len(ids) > 1:
            return {"status": "ambiguous_gene", "query": gene, "candidate_ncbi_ids": ids}
        phenotypes = sorted(self.gene_phenotypes.get(ids[0], {}).values(), key=lambda n: (n["name"], n["index"]))
        return {"status": "ok", "gene_id": ids[0], "namespace": "NCBI",
            "symbols": sorted({n["name"] for n in matches}), "phenotype_count": len(phenotypes),
            "phenotypes": [{"id": n["id"], "namespace": n["source"], "name": n["name"],
                            "node_index": int(n["index"])} for n in phenotypes[:limit]],
            "returned_count": min(limit, len(phenotypes)), "truncated": len(phenotypes) > limit,
            "relation": "phenotype_protein", "evidence_type": "recorded_direct_association",
            "source": "data/processed/primekg/phenotype_disease_gene_edges.csv",
            "model_uses_phenotype_edges": False,
            "warning": "Retrieved graph context, not a model explanation or causal/clinical validation. No recorded association means unknown, not a negative finding."}


def query(repository, method, **kwargs):
    """Keep missing/corrupt local backends distinct from a valid empty result."""
    try:
        return getattr(repository, method)(**kwargs)
    except (DataUnavailable, KeyError, ValueError, TypeError, StopIteration) as error:
        return {"status": "backend_unavailable", "tool": method,
            "message": "A required local result or graph artifact is missing or has an incompatible schema.",
            "detail": str(error) if isinstance(error, DataUnavailable) else "Artifact validation failed"}
