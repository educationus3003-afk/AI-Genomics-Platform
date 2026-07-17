from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests

BASE = Path.home() / "AI_Genomics_Platform"

INPUT_FILE = BASE / "06_Results/Master/database_integrated_table.csv"
OUTPUT_FILE = BASE / "06_Results/Master/database_integrated_table.csv"
CACHE_FILE = BASE / "06_Results/Database_Cache/variant_enrichment_cache.csv"

ENSEMBL_URL = "https://rest.ensembl.org/vep/human/region"
TIMEOUT = 45

OUTPUT_COLUMNS = [
    "dbSNP_ID",
    "gnomAD_AF",
    "gnomAD_Status",
    "AlphaMissense_Score",
    "AlphaMissense_Prediction",
    "SIFT_Prediction",
    "PolyPhen_Prediction",
    "Ensembl_Most_Severe_Consequence",
    "Enrichment_Source",
]


def empty_result(source: str = "Ensembl VEP REST") -> dict[str, Any]:
    return {
        "dbSNP_ID": "Not found",
        "gnomAD_AF": None,
        "gnomAD_Status": "Not evaluated",
        "AlphaMissense_Score": None,
        "AlphaMissense_Prediction": "Not available",
        "SIFT_Prediction": "Not available",
        "PolyPhen_Prediction": "Not available",
        "Ensembl_Most_Severe_Consequence": "Not available",
        "Enrichment_Source": source,
    }


def safe_number(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def query_ensembl(
    chrom: str,
    position: int,
    ref: str,
    alt: str,
) -> dict[str, Any]:
    result = empty_result()

    chrom = str(chrom).replace("chr", "")
    position = int(position)
    ref = str(ref)
    alt = str(alt)

    variant = f"{chrom} {position} . {ref} {alt} . . ."

    payload = {
        "variants": [variant],
    }

    params = {
        "AlphaMissense": "1",
        "dbNSFP": "REVEL_score",
    }

    try:
        response = requests.post(
            ENSEMBL_URL,
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
            },
            params=params,
            json=payload,
            timeout=TIMEOUT,
        )
        response.raise_for_status()
        records = response.json()

        if not records:
            result["Enrichment_Source"] = "Ensembl VEP REST: no record"
            return result

        record = records[0]

        identifiers = record.get("colocated_variants") or []
        rsids: list[str] = []
        frequencies: list[float] = []

        for colocated in identifiers:
            variant_id = colocated.get("id")
            if isinstance(variant_id, str) and variant_id.startswith("rs"):
                rsids.append(variant_id)

            for frequency_group in (colocated.get("frequencies") or {}).values():
                if not isinstance(frequency_group, dict):
                    continue

                for key, value in frequency_group.items():
                    if "gnomad" in key.lower() or key.lower() in {
                        "af",
                        "gnomade",
                        "gnomadg",
                    }:
                        number = safe_number(value)
                        if number is not None:
                            frequencies.append(number)

        if rsids:
            result["dbSNP_ID"] = ";".join(sorted(set(rsids)))

        if frequencies:
            maximum_af = max(frequencies)
            result["gnomAD_AF"] = maximum_af
            result["gnomAD_Status"] = (
                "Common"
                if maximum_af >= 0.05
                else "Rare"
                if maximum_af < 0.01
                else "Low frequency"
            )
        else:
            result["gnomAD_Status"] = "Not found in returned frequencies"

        result["Ensembl_Most_Severe_Consequence"] = record.get(
            "most_severe_consequence",
            "Not available",
        )

        consequences = record.get("transcript_consequences") or []

        canonical = [
            consequence
            for consequence in consequences
            if consequence.get("canonical") == 1
        ]
        selected = canonical or consequences

        alpha_scores: list[float] = []
        alpha_predictions: list[str] = []
        sift_predictions: list[str] = []
        polyphen_predictions: list[str] = []

        for consequence in selected:
            alpha_score = safe_number(
                consequence.get("alphamissense_score")
                or consequence.get("am_score")
            )
            if alpha_score is not None:
                alpha_scores.append(alpha_score)

            alpha_prediction = (
                consequence.get("alphamissense_prediction")
                or consequence.get("am_pathogenicity")
            )
            if alpha_prediction:
                alpha_predictions.append(str(alpha_prediction))

            sift = consequence.get("sift_prediction")
            if sift:
                sift_predictions.append(str(sift))

            polyphen = consequence.get("polyphen_prediction")
            if polyphen:
                polyphen_predictions.append(str(polyphen))

        if alpha_scores:
            score = max(alpha_scores)
            result["AlphaMissense_Score"] = score

            if alpha_predictions:
                result["AlphaMissense_Prediction"] = ";".join(
                    sorted(set(alpha_predictions))
                )
            elif score >= 0.564:
                result["AlphaMissense_Prediction"] = "Likely pathogenic"
            elif score <= 0.34:
                result["AlphaMissense_Prediction"] = "Likely benign"
            else:
                result["AlphaMissense_Prediction"] = "Ambiguous"

        if sift_predictions:
            result["SIFT_Prediction"] = ";".join(
                sorted(set(sift_predictions))
            )

        if polyphen_predictions:
            result["PolyPhen_Prediction"] = ";".join(
                sorted(set(polyphen_predictions))
            )

        return result

    except requests.RequestException as error:
        result["Enrichment_Source"] = f"Ensembl request failed: {error}"
        return result
    except (KeyError, TypeError, ValueError) as error:
        result["Enrichment_Source"] = f"Ensembl parsing failed: {error}"
        return result


def load_cache() -> pd.DataFrame:
    if CACHE_FILE.exists():
        return pd.read_csv(CACHE_FILE)

    return pd.DataFrame(
        columns=["VARIANT_ID", *OUTPUT_COLUMNS]
    )


def main() -> None:
    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Integrated input file does not exist: {INPUT_FILE}"
        )

    dataframe = pd.read_csv(INPUT_FILE)

    required = ["CHROM", "POS", "REF", "ALT"]
    missing = [
        column for column in required
        if column not in dataframe.columns
    ]

    if missing:
        raise ValueError(
            f"Missing required columns: {missing}"
        )

    if "VARIANT_ID" not in dataframe.columns:
        dataframe["VARIANT_ID"] = (
            "chr"
            + dataframe["CHROM"]
            .astype(str)
            .str.replace("chr", "", regex=False)
            + ":"
            + dataframe["POS"].astype(int).astype(str)
            + ":"
            + dataframe["REF"].astype(str)
            + ">"
            + dataframe["ALT"].astype(str)
        )

    cache = load_cache()
    cached_ids = set(cache["VARIANT_ID"].astype(str))

    new_records: list[dict[str, Any]] = []

    print("\nStarting Ensembl/dbSNP/gnomAD/AlphaMissense enrichment...\n")

    for _, row in dataframe.iterrows():
        variant_id = str(row["VARIANT_ID"])

        if variant_id in cached_ids:
            print(f"[CACHE] {variant_id}")
            continue

        print(f"[QUERY] {variant_id}")

        result = query_ensembl(
            chrom=str(row["CHROM"]),
            position=int(row["POS"]),
            ref=str(row["REF"]),
            alt=str(row["ALT"]),
        )

        new_records.append(
            {
                "VARIANT_ID": variant_id,
                **result,
            }
        )

        time.sleep(0.15)

    if new_records:
        cache = pd.concat(
            [cache, pd.DataFrame(new_records)],
            ignore_index=True,
        )

    cache = cache.drop_duplicates(
        subset=["VARIANT_ID"],
        keep="last",
    )

    CACHE_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    cache.to_csv(
        CACHE_FILE,
        index=False,
    )

    dataframe = dataframe.drop(
        columns=[
            column
            for column in OUTPUT_COLUMNS
            if column in dataframe.columns
        ],
        errors="ignore",
    )

    dataframe = dataframe.merge(
        cache[["VARIANT_ID", *OUTPUT_COLUMNS]],
        on="VARIANT_ID",
        how="left",
    )

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    dataframe.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    display_columns = [
        "VARIANT_ID",
        "GENE",
        "dbSNP_ID",
        "gnomAD_AF",
        "AlphaMissense_Score",
        "AlphaMissense_Prediction",
    ]

    display_columns = [
        column
        for column in display_columns
        if column in dataframe.columns
    ]

    print("\nVariant enrichment completed:\n")
    print(
        dataframe[display_columns]
        .to_string(index=False)
    )
    print(f"\nSaved: {OUTPUT_FILE}")
    print(f"Cache: {CACHE_FILE}")


if __name__ == "__main__":
    main()
