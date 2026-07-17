from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests


BASE = Path.home() / "AI_Genomics_Platform"

MASTER_FILE = (
    BASE
    / "06_Results"
    / "Master"
    / "database_integrated_table.csv"
)

INPUT_VCF = (
    BASE
    / "04_Data"
    / "Test"
    / "sample_input.vcf"
)

OUTPUT_FILE = (
    BASE
    / "06_Results"
    / "Master"
    / "database_integrated_table.csv"
)

LIVE_CACHE = (
    BASE
    / "06_Results"
    / "Database_Cache"
    / "live_database_enrichment.json"
)

LIVE_SUMMARY = (
    BASE
    / "06_Results"
    / "Master"
    / "live_database_results.csv"
)

ENSEMBL_URL = (
    "https://rest.ensembl.org/vep/human/region"
)

MYVARIANT_QUERY_URL = (
    "https://myvariant.info/v1/query"
)

TIMEOUT = 50

MISSING_VALUES = {
    "",
    "nan",
    "none",
    "not retrieved",
    "not available",
    "unavailable",
}


def clean(value: Any) -> str:
    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass

    text = str(value).strip()

    if text.lower() in MISSING_VALUES:
        return ""

    return text


def safe_float(value: Any) -> float | None:
    try:
        if value in (None, ""):
            return None

        return float(value)
    except (TypeError, ValueError):
        return None


def parse_info(info_text: str) -> dict[str, str]:
    result: dict[str, str] = {}

    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value

    return result


def read_input_metadata() -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}

    if not INPUT_VCF.exists():
        return records

    with INPUT_VCF.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            chrom, pos, rsid, ref, alt = fields[:5]
            info = parse_info(fields[7])

            chrom = chrom.replace("chr", "")

            variant_id = (
                f"chr{chrom}:{pos}:{ref}>{alt}"
            )

            records[variant_id] = {
                "rsid": rsid if rsid != "." else "",
                "gene": info.get("GENE", ""),
                "consequence": info.get(
                    "CONSEQUENCE",
                    "",
                ),
                "hgvsc": info.get("HGVSC", ""),
                "hgvsp": info.get("HGVSP", ""),
                "clinvar": info.get("CLINVAR", ""),
                "gnomad_af": info.get(
                    "GNOMAD_AF",
                    "",
                ),
                "alphamissense": info.get(
                    "ALPHAMISSENSE",
                    "",
                ),
            }

    return records


def nested_values(
    value: Any,
    key_names: set[str],
) -> list[Any]:
    results: list[Any] = []

    if isinstance(value, dict):
        for key, item in value.items():
            if key.lower() in key_names:
                results.append(item)

            results.extend(
                nested_values(
                    item,
                    key_names,
                )
            )

    elif isinstance(value, list):
        for item in value:
            results.extend(
                nested_values(
                    item,
                    key_names,
                )
            )

    return results


def flatten_text(values: list[Any]) -> list[str]:
    results: list[str] = []

    for value in values:
        if isinstance(value, list):
            results.extend(flatten_text(value))
        elif isinstance(value, dict):
            results.extend(
                flatten_text(list(value.values()))
            )
        else:
            text = clean(value)

            if text:
                results.append(text)

    return results


def query_myvariant(
    *,
    rsid: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> dict[str, Any]:
    result = {
        "clinvar": "",
        "clinvar_review_status": "",
        "dbsnp_id": "",
        "gnomad_af": None,
        "source": "",
        "status": "No MyVariant match",
    }

    queries: list[str] = []

    if rsid.startswith("rs"):
        queries.append(rsid)

    queries.extend(
        [
            f"chr{chrom}:{pos}",
            f"{chrom}:{pos}",
            (
                f"dbsnp.chrom:{chrom} AND "
                f"dbsnp.hg19.start:{pos}"
            ),
        ]
    )

    fields = (
        "clinvar,"
        "dbsnp,"
        "gnomad_exome,"
        "gnomad_genome,"
        "gnomad"
    )

    for query in queries:
        try:
            response = requests.get(
                MYVARIANT_QUERY_URL,
                params={
                    "q": query,
                    "fields": fields,
                    "size": 10,
                },
                timeout=TIMEOUT,
            )

            response.raise_for_status()
            payload = response.json()
            hits = payload.get("hits", [])

            if not hits:
                continue

            selected = None

            for hit in hits:
                hit_text = json.dumps(
                    hit,
                    default=str,
                ).upper()

                if (
                    clean(ref).upper() in hit_text
                    and clean(alt).upper() in hit_text
                ):
                    selected = hit
                    break

            if selected is None:
                selected = hits[0]

            rs_values = flatten_text(
                nested_values(
                    selected.get("dbsnp", {}),
                    {
                        "rsid",
                        "rs",
                        "snp_id",
                    },
                )
            )

            for value in rs_values:
                if value.startswith("rs"):
                    result["dbsnp_id"] = value
                    break

                if value.isdigit():
                    result["dbsnp_id"] = (
                        f"rs{value}"
                    )
                    break

            significance_values = flatten_text(
                nested_values(
                    selected.get("clinvar", {}),
                    {
                        "clinical_significance",
                        "clinsig",
                        "significance",
                    },
                )
            )

            if significance_values:
                result["clinvar"] = "; ".join(
                    sorted(
                        set(significance_values)
                    )
                )

            review_values = flatten_text(
                nested_values(
                    selected.get("clinvar", {}),
                    {
                        "review_status",
                        "review",
                    },
                )
            )

            if review_values:
                result["clinvar_review_status"] = (
                    "; ".join(
                        sorted(set(review_values))
                    )
                )

            frequency_candidates = flatten_text(
                nested_values(
                    {
                        "gnomad_exome": selected.get(
                            "gnomad_exome",
                            {},
                        ),
                        "gnomad_genome": selected.get(
                            "gnomad_genome",
                            {},
                        ),
                        "gnomad": selected.get(
                            "gnomad",
                            {},
                        ),
                    },
                    {
                        "af",
                        "allele_frequency",
                    },
                )
            )

            numeric_frequencies = [
                number
                for number in (
                    safe_float(value)
                    for value in frequency_candidates
                )
                if number is not None
                and 0 <= number <= 1
            ]

            if numeric_frequencies:
                result["gnomad_af"] = max(
                    numeric_frequencies
                )

            result["source"] = "MyVariant.info REST"
            result["status"] = (
                f"MyVariant match using query: {query}"
            )

            return result

        except requests.RequestException:
            continue
        except (
            TypeError,
            ValueError,
            KeyError,
        ):
            continue

    return result


def query_ensembl(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> dict[str, Any]:
    result = {
        "consequence": "",
        "hgvsc": "",
        "hgvsp": "",
        "dbsnp_id": "",
        "alphamissense_score": None,
        "alphamissense_prediction": "",
        "gnomad_af": None,
        "source": "",
        "status": "No Ensembl match",
    }

    variant = (
        f"{chrom} {pos} . {ref} {alt} . . ."
    )

    try:
        response = requests.post(
            ENSEMBL_URL,
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
            },
            params={
                "AlphaMissense": "1",
                "CADD": "1",
                "dbNSFP": "REVEL_score",
                "canonical": "1",
                "hgvs": "1",
            },
            json={
                "variants": [variant],
            },
            timeout=TIMEOUT,
        )

        response.raise_for_status()
        records = response.json()

        if not records:
            return result

        record = records[0]

        result["consequence"] = clean(
            record.get(
                "most_severe_consequence",
                "",
            )
        )

        colocated = (
            record.get("colocated_variants")
            or []
        )

        frequencies: list[float] = []

        for item in colocated:
            identifier = clean(item.get("id"))

            if identifier.startswith("rs"):
                result["dbsnp_id"] = identifier

            frequency_data = (
                item.get("frequencies")
                or {}
            )

            for allele_data in frequency_data.values():
                if not isinstance(
                    allele_data,
                    dict,
                ):
                    continue

                for key, value in allele_data.items():
                    if (
                        "gnomad" in key.lower()
                        or key.lower()
                        in {
                            "af",
                            "gnomade",
                            "gnomadg",
                        }
                    ):
                        number = safe_float(value)

                        if (
                            number is not None
                            and 0 <= number <= 1
                        ):
                            frequencies.append(
                                number
                            )

        if frequencies:
            result["gnomad_af"] = max(
                frequencies
            )

        transcript_consequences = (
            record.get(
                "transcript_consequences"
            )
            or []
        )

        canonical = [
            item
            for item in transcript_consequences
            if item.get("canonical") == 1
        ]

        selected = canonical or transcript_consequences

        alpha_scores: list[float] = []
        alpha_predictions: list[str] = []

        for transcript in selected:
            if not result["hgvsc"]:
                result["hgvsc"] = clean(
                    transcript.get("hgvsc")
                )

            if not result["hgvsp"]:
                result["hgvsp"] = clean(
                    transcript.get("hgvsp")
                )

            score = safe_float(
                transcript.get(
                    "alphamissense_score"
                )
                or transcript.get("am_score")
            )

            if score is not None:
                alpha_scores.append(score)

            prediction = clean(
                transcript.get(
                    "alphamissense_prediction"
                )
                or transcript.get(
                    "am_pathogenicity"
                )
            )

            if prediction:
                alpha_predictions.append(
                    prediction
                )

        if alpha_scores:
            best_score = max(alpha_scores)

            result[
                "alphamissense_score"
            ] = best_score

            if alpha_predictions:
                result[
                    "alphamissense_prediction"
                ] = "; ".join(
                    sorted(
                        set(alpha_predictions)
                    )
                )
            elif best_score >= 0.564:
                result[
                    "alphamissense_prediction"
                ] = "Likely pathogenic"
            elif best_score <= 0.34:
                result[
                    "alphamissense_prediction"
                ] = "Likely benign"
            else:
                result[
                    "alphamissense_prediction"
                ] = "Ambiguous"

        result["source"] = "Ensembl REST VEP"
        result["status"] = "Ensembl REST query successful"

        return result

    except requests.RequestException as error:
        result["status"] = (
            f"Ensembl request failed: {error}"
        )
        return result

    except (
        TypeError,
        ValueError,
        KeyError,
    ) as error:
        result["status"] = (
            f"Ensembl parsing failed: {error}"
        )
        return result


def first_available(
    *values: Any,
    default: str = "Unavailable",
) -> Any:
    for value in values:
        text = clean(value)

        if text:
            return value

    return default


def classify_alpha_applicability(
    consequence: str,
) -> str:
    consequence = clean(
        consequence
    ).lower()

    if "missense" in consequence:
        return "Applicable"

    return (
        "Not applicable — AlphaMissense is "
        "for missense single-nucleotide variants"
    )


def main() -> None:
    if not MASTER_FILE.exists():
        raise FileNotFoundError(
            f"Master table not found: {MASTER_FILE}"
        )

    dataframe = pd.read_csv(MASTER_FILE)
    input_metadata = read_input_metadata()

    cache: dict[str, Any] = {}

    if LIVE_CACHE.exists():
        try:
            cache = json.loads(
                LIVE_CACHE.read_text(
                    encoding="utf-8",
                )
            )
        except json.JSONDecodeError:
            cache = {}

    records: list[dict[str, Any]] = []

    print("\nStarting live database enrichment...\n")

    for _, row in dataframe.iterrows():
        variant_id = clean(
            row.get("VARIANT_ID")
        )

        chrom = clean(
            row.get("CHROM")
        ).replace("chr", "")

        pos = int(row["POS"])
        ref = clean(row.get("REF"))
        alt = clean(row.get("ALT"))
        gene = clean(row.get("GENE"))

        metadata = input_metadata.get(
            variant_id,
            {},
        )

        rsid = first_available(
            row.get("dbSNP_ID"),
            row.get("Final_dbSNP_ID"),
            metadata.get("rsid"),
            default="",
        )

        if variant_id in cache:
            combined = cache[variant_id]
            print(f"[CACHE] {variant_id}")
        else:
            print(f"[LIVE] {variant_id}")

            ensembl = query_ensembl(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
            )

            myvariant = query_myvariant(
                rsid=clean(rsid),
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
            )

            consequence = first_available(
                ensembl.get("consequence"),
                row.get("VEP_Consequence"),
                metadata.get("consequence"),
            )

            alpha_applicability = (
                classify_alpha_applicability(
                    consequence
                )
            )

            if (
                alpha_applicability
                != "Applicable"
            ):
                alpha_prediction = (
                    alpha_applicability
                )
                alpha_score = None
            else:
                alpha_prediction = (
                    first_available(
                        ensembl.get(
                            "alphamissense_prediction"
                        ),
                        metadata.get(
                            "alphamissense"
                        ),
                    )
                )

                alpha_score = ensembl.get(
                    "alphamissense_score"
                )

            combined = {
                "VARIANT_ID": variant_id,
                "GENE": gene,
                "Live_dbSNP_ID": first_available(
                    ensembl.get("dbsnp_id"),
                    myvariant.get("dbsnp_id"),
                    rsid,
                ),
                "Live_ClinVar": first_available(
                    myvariant.get("clinvar"),
                    metadata.get("clinvar"),
                ),
                "Live_ClinVar_ReviewStatus": (
                    first_available(
                        myvariant.get(
                            "clinvar_review_status"
                        )
                    )
                ),
                "Live_gnomAD_AF": first_available(
                    ensembl.get("gnomad_af"),
                    myvariant.get("gnomad_af"),
                    metadata.get("gnomad_af"),
                ),
                "Live_VEP_Consequence": (
                    consequence
                ),
                "Live_VEP_HGVSc": first_available(
                    ensembl.get("hgvsc"),
                    row.get("VEP_HGVSc"),
                    metadata.get("hgvsc"),
                ),
                "Live_VEP_HGVSp": first_available(
                    ensembl.get("hgvsp"),
                    row.get("VEP_HGVSp"),
                    metadata.get("hgvsp"),
                ),
                "Live_AlphaMissense_Score": (
                    alpha_score
                ),
                "Live_AlphaMissense_Prediction": (
                    alpha_prediction
                ),
                "AlphaMissense_Applicability": (
                    alpha_applicability
                ),
                "Ensembl_Status": ensembl.get(
                    "status"
                ),
                "MyVariant_Status": myvariant.get(
                    "status"
                ),
                "Evidence_Sources": "; ".join(
                    source
                    for source in [
                        clean(
                            ensembl.get("source")
                        ),
                        clean(
                            myvariant.get("source")
                        ),
                        (
                            "Validated input VCF metadata"
                            if metadata
                            else ""
                        ),
                    ]
                    if source
                ),
            }

            cache[variant_id] = combined

            time.sleep(0.35)

        records.append(combined)

    LIVE_CACHE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    LIVE_CACHE.write_text(
        json.dumps(
            cache,
            indent=2,
            default=str,
        ),
        encoding="utf-8",
    )

    live = pd.DataFrame(records)

    live.to_csv(
        LIVE_SUMMARY,
        index=False,
    )

    live_columns = [
        column
        for column in live.columns
        if column
        not in {
            "GENE",
        }
    ]

    dataframe = dataframe.drop(
        columns=[
            column
            for column in live_columns
            if column in dataframe.columns
            and column != "VARIANT_ID"
        ],
        errors="ignore",
    )

    dataframe = dataframe.merge(
        live[live_columns],
        on="VARIANT_ID",
        how="left",
    )

    dataframe.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    display = [
        "VARIANT_ID",
        "GENE",
        "Live_dbSNP_ID",
        "Live_ClinVar",
        "Live_gnomAD_AF",
        "Live_VEP_Consequence",
        "Live_AlphaMissense_Prediction",
        "Evidence_Sources",
    ]

    display = [
        column
        for column in display
        if column in dataframe.columns
    ]

    print("\nLive enrichment completed:\n")

    print(
        dataframe[display]
        .to_string(index=False)
    )

    print(f"\nSaved: {OUTPUT_FILE}")
    print(f"Live summary: {LIVE_SUMMARY}")
    print(f"Cache: {LIVE_CACHE}")


if __name__ == "__main__":
    main()
