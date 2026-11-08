"""Local catalog only; no KG or prediction backend is connected yet."""

from agents import function_tool

SUBTYPES = ("GI.CIN", "GI.HM-indel", "GI.HM-SNV", "GI.GS")


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
            "source": "Local ECoGraph subtype catalog; no GNN predictions are available"}
