import json
from pathlib import Path

import pandas as pd


BASE = Path.home() / "AI_Genomics_Platform"

INPUT_FILE = (
    BASE / "06_Results/Master/database_integrated_table.csv"
)

OUTPUT_FILE = (
    BASE / "06_Results/Master/integrated_evidence.json"
)


def clean(value):
    if pd.isna(value):
        return "Not available"

    text = str(value).strip()

    if text in {
        "",
        "nan",
        "NaN",
        "Not Retrieved",
        "Not available",
    }:
        return "Not available"

    return text


def main():
    dataframe = pd.read_csv(INPUT_FILE)

    records = []

    for _, row in dataframe.iterrows():
        record = {
            "variant_id": clean(row.get("VARIANT_ID")),
            "gene": clean(row.get("GENE")),
            "coordinates": {
                "chromosome": clean(row.get("CHROM")),
                "position": clean(row.get("POS")),
                "reference": clean(row.get("REF")),
                "alternate": clean(row.get("ALT")),
            },
            "intervar": {
                "classification": clean(
                    row.get("INTERVAR_CLASSIFICATION")
                ),
                "acmg_evidence": clean(
                    row.get("ACMG_EVIDENCE")
                ),
                "source": clean(row.get("SOURCE")),
            },
            "clinvar": {
                "clinical_significance": clean(
                    row.get("ClinVar_ClinicalSignificance")
                ),
                "review_status": clean(
                    row.get("ClinVar_ReviewStatus")
                ),
                "disease": clean(
                    row.get("ClinVar_Disease")
                ),
                "variation_id": clean(
                    row.get("ClinVar_VariationID")
                ),
            },
            "gnomad": {
                "allele_frequency": clean(
                    row.get("gnomAD_AlleleFrequency")
                ),
                "allele_count": clean(
                    row.get("gnomAD_AlleleCount")
                ),
                "population": clean(
                    row.get("gnomAD_Population")
                ),
            },
            "predictions": {
                "alphamissense_score": clean(
                    row.get("AlphaMissense_Score")
                ),
                "alphamissense_prediction": clean(
                    row.get("AlphaMissense_Prediction")
                ),
                "revel_score": clean(
                    row.get("REVEL_Score")
                ),
                "cadd_raw": clean(
                    row.get("CADD_RawScore")
                ),
                "cadd_phred": clean(
                    row.get("CADD_PHRED")
                ),
                "spliceai_delta_score": clean(
                    row.get("SpliceAI_DeltaScore")
                ),
            },
            "disease_resources": {
                "omim_disease": clean(
                    row.get("OMIM_DiseaseName")
                ),
                "omim_inheritance": clean(
                    row.get("OMIM_Inheritance")
                ),
                "omim_id": clean(
                    row.get("OMIM_OMIM_ID")
                ),
                "uniprot_protein": clean(
                    row.get("UniProt_Protein")
                ),
                "uniprot_function": clean(
                    row.get("UniProt_Function")
                ),
                "clingen_validity": clean(
                    row.get("ClinGen_DiseaseValidity")
                ),
            },
            "vep": {
                "consequence": clean(
                    row.get("VEP_Consequence")
                ),
                "hgvsc": clean(row.get("VEP_HGVSc")),
                "hgvsp": clean(row.get("VEP_HGVSp")),
            },
        }

        records.append(record)

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_FILE.open("w", encoding="utf-8") as file:
        json.dump(records, file, indent=2)

    print(f"Exported {len(records)} variant records.")
    print(f"Saved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
