import json
from pathlib import Path
from typing import Any


BASE = Path.home() / "AI_Genomics_Platform"
KNOWLEDGE_FILE = BASE / "01_Disease_Knowledge/diseases.json"


def load_knowledge_base() -> dict[str, Any]:
    if not KNOWLEDGE_FILE.exists():
        raise FileNotFoundError(f"Knowledge base not found: {KNOWLEDGE_FILE}")

    with KNOWLEDGE_FILE.open("r", encoding="utf-8") as file:
        return json.load(file)


def get_disease_profile(gene: str) -> dict[str, Any]:
    knowledge = load_knowledge_base()
    gene = gene.strip().upper()

    if gene not in knowledge:
        return {
            "gene": gene,
            "disease": "Unknown",
            "status": "No disease profile currently available"
        }

    return knowledge[gene]


def main() -> None:
    knowledge = load_knowledge_base()

    print("\nAvailable disease profiles:\n")

    for gene, profile in knowledge.items():
        print(f"{gene}: {profile['disease']}")
        print(f"  Inheritance: {profile['inheritance']}")
        print(f"  Mechanism: {profile['mechanism']}")
        print()


if __name__ == "__main__":
    main()
