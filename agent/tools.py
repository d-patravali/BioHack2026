"""Tools for local graph evidence and explicitly cached GNN experiment outputs."""

from agents import function_tool
from functools import lru_cache

from .repository import ResearchRepository, SUBTYPES, query


@lru_cache(maxsize=1)
def repository():
    return ResearchRepository()


@function_tool
def get_subtypes(disease: str) -> dict:
    """Retrieve the supported disease-subtype catalog.

    Args:
        disease: Disease identifier. Currently only colorectal_cancer is supported.
    """
    normalized = disease.strip().lower().replace(" ", "_").replace("-", "_")
    if normalized != "colorectal_cancer":
        return {"status": "unsupported_disease", "disease": normalized, "subtypes": [],
                "supported_diseases": ["colorectal_cancer"]}
    return {"status": "ok", "disease": normalized, "subtypes": list(SUBTYPES),
            "source": "Local ECoGraph subtype catalog; cached GNN experiments available, no fresh inference"}


@function_tool
def get_drug_predictions(subtype: str, top_k: int = 10) -> dict:
    """Read cached experimental GNN scores for a CRC subtype, including all cutoff ties.

    Args:
        subtype: GI.CIN, GI.HM-indel, GI.HM-SNV, or GI.GS.
        top_k: Requested cutoff from 1 to 20; full tied groups can exceed this count.
    """
    return query(repository(), "get_drug_predictions", subtype=subtype, top_k=top_k)


@function_tool
def get_drug_targets(drug_id: str) -> dict:
    """Retrieve recorded direct PrimeKG drug targets, excluding other drug-protein roles.

    Args:
        drug_id: Exact DrugBank ID such as DB00596.
    """
    return query(repository(), "get_drug_targets", drug_id=drug_id)


@function_tool
def get_subtype_genes(subtype: str, limit: int = 20, gene: str = "") -> dict:
    """Retrieve the most frequent observed subtype genes, with cohort counts and frequencies.

    Args:
        subtype: One of the four supported CRC subtype identifiers.
        limit: Maximum genes returned, from 1 to 50. Remaining genes are explicitly truncated.
        gene: Optional exact symbol or NCBI ID to check across the full profile, not just the displayed top genes.
    """
    return query(repository(), "get_subtype_genes", subtype=subtype, limit=limit, gene=gene)


@function_tool
def get_gene_phenotypes(gene: str, limit: int = 10) -> dict:
    """Query explicit PrimeKG gene-phenotype associations as context, not model explanations.

    Args:
        gene: NCBI gene identifier or exact recorded symbol, e.g. 2908 or NR3C1.
        limit: Maximum associations returned, from 1 to 50.
    """
    return query(repository(), "get_gene_phenotypes", gene=gene, limit=limit)


@function_tool
def get_model_evaluation() -> dict:
    """Retrieve held-out AP, popularity comparison, and molecular-context ablation results.

    These evaluate sampled disease-indication recovery, not subtype therapeutic efficacy.
    """
    return query(repository(), "get_model_evaluation")


ORIGINAL_TOOLS = (get_subtypes, get_drug_predictions, get_drug_targets,
                  get_subtype_genes, get_gene_phenotypes, get_model_evaluation)


def tools_for_backend(backend='ecograph'):
    """Bind each agent's tools to a fixed backend, preserving original defaults."""
    from .backends import make_repository, validate_backend
    validate_backend(backend)
    if backend == 'ecograph':
        return list(ORIGINAL_TOOLS)
    selected = make_repository(backend)

    @function_tool
    def get_subtypes(disease: str) -> dict:
        """List CRCSC expression-defined CMS1–CMS4 for colorectal_cancer."""
        return query(selected, 'get_subtypes', disease=disease)

    @function_tool
    def get_drug_predictions(subtype: str, top_k: int = 10) -> dict:
        """Read cached CRCSC signed-profile CMS1–4 hypotheses, retaining all cutoff ties. top_k: 1–20."""
        return query(selected, 'get_drug_predictions', subtype=subtype, top_k=top_k)

    @function_tool
    def get_drug_targets(drug_id: str) -> dict:
        """Retrieve recorded direct PrimeKG targets for a DrugBank ID, excluding other roles."""
        return query(selected, 'get_drug_targets', drug_id=drug_id)

    @function_tool
    def get_subtype_genes(subtype: str, limit: int = 20, gene: str = '') -> dict:
        """Retrieve signed CMS expression associations, NOT mutation counts. limit: 1–50; gene: optional exact symbol/NCBI ID."""
        return query(selected, 'get_subtype_genes', subtype=subtype, limit=limit, gene=gene)

    @function_tool
    def get_gene_phenotypes(gene: str, limit: int = 10) -> dict:
        """Retrieve explicit PrimeKG phenotype–gene context, not model explanations. limit: 1–50."""
        return query(selected, 'get_gene_phenotypes', gene=gene, limit=limit)

    @function_tool
    def get_model_evaluation() -> dict:
        """Read the CRCSC experiment's actual popularity, PPI/pathway and drug-only AP; not subtype efficacy."""
        return query(selected, 'get_model_evaluation')

    return [get_subtypes, get_drug_predictions, get_drug_targets,
            get_subtype_genes, get_gene_phenotypes, get_model_evaluation]
