from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


BASE = Path.home() / "AI_Genomics_Platform"
ANNOTATION_FILE = BASE / "06_Results/Master/standardized_annotation_table.csv"
KNOWLEDGE_FILE = BASE / "01_Disease_Knowledge/diseases.json"
OUTPUT_JSON = BASE / "06_Results/Interpretation/ai_interpretations.json"
OUTPUT_CSV = BASE / "06_Results/Interpretation/ai_interpretations.csv"


def load_annotations() -> pd.DataFrame:
    if not ANNOTATION_FILE.exists():
        raise FileNotFoundError(f"Annotation file not found: {ANNOTATION_FILE}")

    return pd.read_csv(ANNOTATION_FILE)


def load_knowledge() -> dict[str, Any]:
    if not KNOWLEDGE_FILE.exists():
        raise FileNotFoundError(f"Knowledge file not found: {KNOWLEDGE_FILE}")

    with KNOWLEDGE_FILE.open("r", encoding="utf-8") as file:
        return json.load(file)


def clean_text(value: Any, default: str = "Not available") -> str:
    if pd.isna(value):
        return default

    text = str(value).strip()
    return text if text else default


def infer_variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "single nucleotide variant"

    if len(ref) > len(alt):
        return "deletion"

    if len(alt) > len(ref):
        return "insertion"

    return "complex variant"


def build_interpretation(
    row: pd.Series,
    knowledge: dict[str, Any],
) -> dict[str, Any]:
    gene = clean_text(row.get("GENE")).upper()
    profile = knowledge.get(gene, {})

    disease = profile.get("disease", "Disease association not available")
    inheritance = profile.get("inheritance", "Not available")
    mechanism = profile.get("mechanism", "Not available")
    protein = profile.get("protein", "Not available")
    features = profile.get("major_features", [])
    validation_focus = profile.get("validation_focus", [])

    variant_id = clean_text(row.get("VARIANT_ID"))
    ref = clean_text(row.get("REF"))
    alt = clean_text(row.get("ALT"))
    variant_type = infer_variant_type(ref, alt)

    intervar_class = clean_text(
        row.get("INTERVAR_CLASSIFICATION"),
        "Pending expert review",
    )

    acmg_evidence = clean_text(
        row.get("ACMG_EVIDENCE"),
        "Pending expert review",
    )

    source = clean_text(row.get("SOURCE"))

    gene_summary = (
        f"{gene} encodes {protein}. "
        f"The primary disease mechanism associated with this case study is: "
        f"{mechanism}."
    )

    variant_summary = (
        f"The detected variant {variant_id} is a {variant_type} affecting "
        f"{gene}. The current automated classification is "
        f"{intervar_class}."
    )

    disease_summary = (
        f"{gene} is associated with {disease}, which is typically inherited "
        f"in an {inheritance.lower()} pattern."
    )

    clinical_features = (
        ", ".join(features)
        if features
        else "Clinical features are not available in the current knowledge base."
    )

    clinical_interpretation = (
        f"This variant should be interpreted together with the patient's "
        f"phenotype, inheritance pattern, population frequency, ClinVar status, "
        f"functional prediction scores, transcript context, and ACMG evidence. "
        f"Relevant clinical features for {disease} include: {clinical_features}"
    )

    acmg_summary = (
        f"Automated ACMG evidence: {acmg_evidence}. "
        f"This result is decision support only and requires review by a "
        f"qualified genetics professional."
    )

    patient_explanation = (
        f"A change was detected in the {gene} gene. This gene is linked to "
        f"{disease}. The result does not by itself confirm a diagnosis and "
        f"must be considered together with symptoms, family history, and "
        f"clinical testing."
    )

    return {
        "variant_id": variant_id,
        "gene": gene,
        "disease": disease,
        "inheritance": inheritance,
        "variant_type": variant_type,
        "automated_classification": intervar_class,
        "acmg_evidence": acmg_evidence,
        "annotation_source": source,
        "gene_summary": gene_summary,
        "variant_summary": variant_summary,
        "disease_summary": disease_summary,
        "clinical_interpretation": clinical_interpretation,
        "acmg_summary": acmg_summary,
        "patient_friendly_explanation": patient_explanation,
        "validation_focus": validation_focus,
        "clinical_review_required": True,
    }


def main() -> None:
    annotations = load_annotations()
    knowledge = load_knowledge()

    interpretations = [
        build_interpretation(row, knowledge)
        for _, row in annotations.iterrows()
    ]

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_JSON.open("w", encoding="utf-8") as file:
        json.dump(interpretations, file, indent=2, ensure_ascii=False)

    pd.DataFrame(interpretations).to_csv(OUTPUT_CSV, index=False)

    print(f"\nGenerated {len(interpretations)} interpretations.")
    print(f"JSON: {OUTPUT_JSON}")
    print(f"CSV:  {OUTPUT_CSV}")

    for item in interpretations:
        print("\n" + "=" * 80)
        print(f"Variant: {item['variant_id']}")
        print(f"Gene: {item['gene']}")
        print(f"Disease: {item['disease']}")
        print(f"Classification: {item['automated_classification']}")
        print("\nClinical interpretation:")
        print(item["clinical_interpretation"])


if __name__ == "__main__":
    main()
