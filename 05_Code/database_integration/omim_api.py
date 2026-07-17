"""
GeneInsightAI - OMIM Integration

This module provides demonstration OMIM metadata for the
project's example diseases. If an OMIM API key is configured
in the future, this module can be extended to perform live
queries.
"""

from typing import Dict, Optional

OMIM_DATA: Dict[str, Dict[str, str]] = {
    "CFTR": {
        "gene": "CFTR",
        "disease": "Cystic fibrosis",
        "omim_id": "219700",
        "inheritance": "Autosomal recessive",
        "status": "Reference metadata",
    },
    "HBB": {
        "gene": "HBB",
        "disease": "Sickle cell disease",
        "omim_id": "603903",
        "inheritance": "Autosomal recessive",
        "status": "Reference metadata",
    },
    "FBN1": {
        "gene": "FBN1",
        "disease": "Marfan syndrome",
        "omim_id": "154700",
        "inheritance": "Autosomal dominant",
        "status": "Reference metadata",
    },
    "MECP2": {
        "gene": "MECP2",
        "disease": "Rett syndrome",
        "omim_id": "312750",
        "inheritance": "X-linked dominant",
        "status": "Reference metadata",
    },
}


def get_omim_information(gene: str) -> Dict[str, Optional[str]]:
    """
    Return OMIM reference information for a gene.
    """

    gene = gene.upper()

    if gene in OMIM_DATA:
        return OMIM_DATA[gene]

    return {
        "gene": gene,
        "disease": "Not available",
        "omim_id": "Not available",
        "inheritance": "Not available",
        "status": (
            "Gene not found in local reference. "
            "Live OMIM API integration not configured."
        ),
    }


if __name__ == "__main__":
    for gene in ["CFTR", "HBB", "FBN1", "MECP2", "BRCA1"]:
        print(get_omim_information(gene))
