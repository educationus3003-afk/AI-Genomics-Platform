from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests


BASE = Path.home() / "AI_Genomics_Platform"

INPUT_FILE = BASE / "06_Results/Master/database_integrated_table.csv"
OUTPUT_FILE = BASE / "06_Results/Master/database_integrated_table.csv"
CACHE_FILE = BASE / "04_Data/Database_Resources/UniProt/uniprot_results.csv"

API_URL = "https://rest.uniprot.org/uniprotkb/search"

TARGET_GENES = {"CFTR", "HBB", "MECP2", "FBN1"}


def clean(value: Any, default: str = "Not available") -> str:
    if value is None:
        return default

    text = str(value).strip()

    if text in {"", ".", "nan", "NaN", "None"}:
        return default

    return text


def get_recommended_name(entry: dict[str, Any]) -> str:
    description = entry.get("proteinDescription", {})

    recommended = description.get("recommendedName", {})
    full_name = recommended.get("fullName", {})
    value = full_name.get("value")

    if value:
        return clean(value)

    submission = description.get("submissionNames", [])

    if submission:
        return clean(
            submission[0].get("fullName", {}).get("value")
        )

    return "Not available"


def get_function(entry: dict[str, Any]) -> str:
    comments = entry.get("comments", [])

    function_texts: list[str] = []

    for comment in comments:
        if comment.get("commentType") != "FUNCTION":
            continue

        for text_block in comment.get("texts", []):
            value = text_block.get("value")

            if value:
                function_texts.append(clean(value))

    if not function_texts:
        return "Not available"

    return " ".join(function_texts)


def get_domains(entry: dict[str, Any]) -> str:
    features = entry.get("features", [])

    domain_names: list[str] = []

    for feature in features:
        feature_type = feature.get("type", "")

        if feature_type not in {
            "Domain",
            "Region",
            "Repeat",
            "Motif",
        }:
            continue

        description = feature.get("description")

        if description:
            domain_names.append(clean(description))

    unique_domains = list(dict.fromkeys(domain_names))

    if not unique_domains:
        return "Not available"

    return "; ".join(unique_domains[:20])


def query_uniprot(gene: str) -> dict[str, str]:
    params = {
        "query": (
            f"gene_exact:{gene} AND "
            "organism_id:9606 AND reviewed:true"
        ),
        "format": "json",
        "size": 1,
    }

    response = requests.get(
        API_URL,
        params=params,
        timeout=60,
        headers={
            "User-Agent": "AI-Genomics-Platform/1.0"
        },
    )

    response.raise_for_status()

    payload = response.json()
    results = payload.get("results", [])

    if not results:
        return {
            "GENE": gene,
            "UniProt_Accession": "Not available",
            "UniProt_Protein": "Not available",
            "UniProt_Function": "Not available",
            "UniProt_Domains": "Not available",
            "UniProt_SequenceLength": "Not available",
            "UniProt_Source": "UniProt REST API",
        }

    entry = results[0]

    sequence_length = (
        entry.get("sequence", {}).get("length")
    )

    return {
        "GENE": gene,
        "UniProt_Accession": clean(
            entry.get("primaryAccession")
        ),
        "UniProt_Protein": get_recommended_name(entry),
        "UniProt_Function": get_function(entry),
        "UniProt_Domains": get_domains(entry),
        "UniProt_SequenceLength": clean(sequence_length),
        "UniProt_Source": "UniProt REST API",
    }


def main() -> None:
    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Integrated table not found: {INPUT_FILE}"
        )

    dataframe = pd.read_csv(INPUT_FILE)

    genes = sorted(
        {
            str(gene).strip().upper()
            for gene in dataframe["GENE"].dropna()
            if str(gene).strip().upper() in TARGET_GENES
        }
    )

    records: list[dict[str, str]] = []

    for gene in genes:
        print(f"Querying UniProt for {gene}...")

        try:
            records.append(query_uniprot(gene))
        except requests.RequestException as error:
            records.append(
                {
                    "GENE": gene,
                    "UniProt_Accession": "API error",
                    "UniProt_Protein": "API error",
                    "UniProt_Function": f"API error: {error}",
                    "UniProt_Domains": "API error",
                    "UniProt_SequenceLength": "API error",
                    "UniProt_Source": "UniProt REST API",
                }
            )

        time.sleep(0.3)

    uniprot = pd.DataFrame(records)

    CACHE_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    uniprot.to_csv(
        CACHE_FILE,
        index=False,
    )

    columns_to_replace = [
        "UniProt_Accession",
        "UniProt_Protein",
        "UniProt_Function",
        "UniProt_Domains",
        "UniProt_SequenceLength",
        "UniProt_Source",
    ]

    dataframe = dataframe.drop(
        columns=[
            column
            for column in columns_to_replace
            if column in dataframe.columns
        ],
        errors="ignore",
    )

    dataframe = dataframe.merge(
        uniprot,
        on="GENE",
        how="left",
    )

    dataframe.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    display_columns = [
        "GENE",
        "UniProt_Accession",
        "UniProt_Protein",
        "UniProt_SequenceLength",
        "UniProt_Function",
    ]

    print("\nUniProt integration complete:\n")
    print(
        dataframe[display_columns]
        .drop_duplicates(subset=["GENE"])
        .to_string(index=False)
    )

    print(f"\nSaved: {OUTPUT_FILE}")
    print(f"Cache: {CACHE_FILE}")


if __name__ == "__main__":
    main()
