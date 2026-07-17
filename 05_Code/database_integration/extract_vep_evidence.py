from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


BASE = Path.home() / "AI_Genomics_Platform"

VEP_VCF = BASE / "06_Results/VEP/all_four_diseases_vep.vcf"
MASTER_TABLE = BASE / "06_Results/Master/database_integrated_table.csv"
KNOWLEDGE_FILE = BASE / "01_Disease_Knowledge/diseases.json"

OUTPUT_FILE = BASE / "06_Results/Master/database_integrated_table.csv"
EXTRACTED_FILE = BASE / "06_Results/tables/vep_extracted_evidence.csv"


def parse_info(info_text: str) -> dict[str, str]:
    info: dict[str, str] = {}

    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        elif item:
            info[item] = "True"

    return info


def clean_value(value: Any, default: str = "Not Retrieved") -> str:
    if value is None:
        return default

    text = str(value).strip()

    if text in {"", ".", "nan", "NaN", "None"}:
        return default

    return text


def review_status_from_stars(value: str) -> str:
    value = clean_value(value)

    if value == "Not Retrieved":
        return value

    mapping = {
        "0": "No assertion criteria provided",
        "1": "Criteria provided, single submitter",
        "2": "Criteria provided, multiple submitters",
        "3": "Reviewed by expert panel",
        "4": "Practice guideline",
    }

    return mapping.get(value, f"{value} ClinVar review stars")


def alpha_prediction(value: str) -> str:
    value = clean_value(value)

    if value == "Not Retrieved":
        return value

    lowered = value.lower()

    if "pathogenic" in lowered or "benign" in lowered or "ambiguous" in lowered:
        return value

    try:
        score = float(value)

        if score >= 0.564:
            return "Likely pathogenic"
        if score <= 0.34:
            return "Likely benign"
        return "Ambiguous"
    except ValueError:
        return value


def load_knowledge() -> dict[str, Any]:
    if not KNOWLEDGE_FILE.exists():
        return {}

    with KNOWLEDGE_FILE.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def extract_vep_records() -> pd.DataFrame:
    records: list[dict[str, Any]] = []

    if not VEP_VCF.exists():
        raise FileNotFoundError(f"VEP VCF not found: {VEP_VCF}")

    knowledge = load_knowledge()

    with VEP_VCF.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            chrom, pos, variant_id, ref, alt, qual, filt, info_text = fields[:8]
            info = parse_info(info_text)

            gene = clean_value(info.get("GENE"), "Unknown").upper()

            profile = knowledge.get(gene, {})
            protein = clean_value(profile.get("protein"))

            for alternate in alt.split(","):
                normalized_chrom = chrom.replace("chr", "")

                records.append(
                    {
                        "VARIANT_ID": (
                            f"chr{normalized_chrom}:{pos}:{ref}>{alternate}"
                        ),
                        "CHROM": normalized_chrom,
                        "POS": int(pos),
                        "REF": ref,
                        "ALT": alternate,
                        "GENE": gene,

                        "ClinVar_ClinicalSignificance":
                            clean_value(info.get("CLINVAR")),

                        "ClinVar_ReviewStatus":
                            review_status_from_stars(
                                info.get("CLINVAR_STARS", "")
                            ),

                        "ClinVar_Disease":
                            clean_value(
                                info.get("CLNDN", info.get("OMIM"))
                            ),

                        "ClinVar_VariationID":
                            clean_value(
                                info.get("CLINVAR_VARIATION_ID")
                            ),

                        "gnomAD_AlleleFrequency":
                            clean_value(info.get("GNOMAD_AF")),

                        "gnomAD_AlleleCount":
                            clean_value(info.get("GNOMAD_AC")),

                        "gnomAD_Population":
                            clean_value(info.get("GNOMAD_POPMAX")),

                        "OMIM_DiseaseName":
                            clean_value(info.get("OMIM")),

                        "OMIM_Inheritance":
                            clean_value(info.get("INHERITANCE")),

                        "OMIM_OMIM_ID":
                            clean_value(info.get("OMIM_ID")),

                        "UniProt_Protein": protein,

                        "AlphaMissense_Score":
                            clean_value(info.get("ALPHAMISSENSE")),

                        "AlphaMissense_Prediction":
                            alpha_prediction(
                                info.get("ALPHAMISSENSE", "")
                            ),

                        "REVEL_Score":
                            clean_value(info.get("REVEL")),

                        "CADD_RawScore":
                            clean_value(info.get("CADD_RAW")),

                        "CADD_PHRED":
                            clean_value(info.get("CADD_PHRED")),

                        "SpliceAI_DeltaScore":
                            clean_value(
                                info.get(
                                    "SPLICEAI_MAX",
                                    info.get("SpliceAI")
                                )
                            ),

                        "VEP_Consequence":
                            clean_value(info.get("CONSEQUENCE")),

                        "VEP_HGVSc":
                            clean_value(info.get("HGVSC")),

                        "VEP_HGVSp":
                            clean_value(info.get("HGVSP")),

                        "VEP_Source":
                            "VEP release 116 offline annotation",
                    }
                )

    return pd.DataFrame(records)


def main() -> None:
    if not MASTER_TABLE.exists():
        raise FileNotFoundError(
            f"Database table not found: {MASTER_TABLE}"
        )

    master = pd.read_csv(MASTER_TABLE)
    evidence = extract_vep_records()

    EXTRACTED_FILE.parent.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(EXTRACTED_FILE, index=False)

    keys = ["CHROM", "POS", "REF", "ALT", "GENE"]

    master["CHROM"] = master["CHROM"].astype(str).str.replace(
        "chr", "", regex=False
    )
    evidence["CHROM"] = evidence["CHROM"].astype(str)

    master["POS"] = pd.to_numeric(master["POS"], errors="coerce")
    evidence["POS"] = pd.to_numeric(evidence["POS"], errors="coerce")

    evidence = evidence.drop_duplicates(subset=keys, keep="first")

    result = master.copy()

    evidence_columns = [
        column for column in evidence.columns
        if column not in {"VARIANT_ID", *keys}
    ]

    merged = master.merge(
        evidence[keys + evidence_columns],
        on=keys,
        how="left",
        suffixes=("", "_VEP"),
    )

    # Permit mixed text and numeric evidence during column updates.
    # This avoids pandas ArrowStringArray/int64 assignment errors.
    merged = merged.astype(object)

    for column in evidence_columns:
        vep_column = f"{column}_VEP"

        if vep_column in merged.columns:
            existing = merged[column].astype(str)

            replace_mask = (
                merged[column].isna()
                | existing.isin(
                    [
                        "Not Retrieved",
                        "Not available",
                        "nan",
                        "NaN",
                        "",
                    ]
                )
            )

            merged.loc[replace_mask, column] = merged.loc[
                replace_mask, vep_column
            ]

            merged = merged.drop(columns=[vep_column])

    merged.to_csv(OUTPUT_FILE, index=False)

    print("\nReal VEP evidence integrated successfully.\n")

    display_columns = [
        "VARIANT_ID",
        "GENE",
        "ClinVar_ClinicalSignificance",
        "ClinVar_ReviewStatus",
        "gnomAD_AlleleFrequency",
        "AlphaMissense_Prediction",
        "REVEL_Score",
        "CADD_PHRED",
        "SpliceAI_DeltaScore",
        "OMIM_DiseaseName",
        "UniProt_Protein",
    ]

    available = [
        column for column in display_columns
        if column in merged.columns
    ]

    print(merged[available].to_string(index=False))
    print(f"\nSaved: {OUTPUT_FILE}")
    print(f"Extracted evidence: {EXTRACTED_FILE}")


if __name__ == "__main__":
    main()
