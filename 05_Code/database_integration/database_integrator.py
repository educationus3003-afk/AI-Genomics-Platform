from pathlib import Path
import json

import pandas as pd

from omim_api import get_omim_information


BASE = Path.home() / "AI_Genomics_Platform"

ANNOTATION_FILE = (
    BASE / "06_Results/Master/standardized_annotation_table.csv"
)
DATABASE_SCHEMA = (
    BASE / "04_Data/Database_Resources/database_schema.json"
)
OUTPUT_FILE = (
    BASE / "06_Results/Master/database_integrated_table.csv"
)


def load_schema() -> dict:
    if not DATABASE_SCHEMA.exists():
        raise FileNotFoundError(
            f"Schema not found: {DATABASE_SCHEMA}"
        )

    with DATABASE_SCHEMA.open("r", encoding="utf-8") as file:
        return json.load(file)


def add_omim_columns(dataframe: pd.DataFrame) -> pd.DataFrame:
    result = dataframe.copy()

    if "GENE" not in result.columns:
        result["OMIM_OMIM_ID"] = "Not available"
        result["OMIM_DiseaseName"] = "Not available"
        result["OMIM_Inheritance"] = "Not available"
        result["OMIM_Status"] = "GENE column missing"
        return result

    omim_records = result["GENE"].fillna("").astype(str).apply(
        get_omim_information
    )

    result["OMIM_OMIM_ID"] = omim_records.apply(
        lambda item: item["omim_id"]
    )
    result["OMIM_DiseaseName"] = omim_records.apply(
        lambda item: item["disease"]
    )
    result["OMIM_Inheritance"] = omim_records.apply(
        lambda item: item["inheritance"]
    )
    result["OMIM_Status"] = omim_records.apply(
        lambda item: item["status"]
    )

    return result


def main() -> None:
    if not ANNOTATION_FILE.exists():
        raise FileNotFoundError(
            f"Annotation file not found: {ANNOTATION_FILE}"
        )

    dataframe = pd.read_csv(ANNOTATION_FILE)
    schema = load_schema()

    for database, fields in schema.items():
        for field in fields:
            column_name = f"{database}_{field}"

            if column_name not in dataframe.columns:
                dataframe[column_name] = "Not Retrieved"

    dataframe = add_omim_columns(dataframe)

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    dataframe.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    display_columns = [
        "GENE",
        "OMIM_OMIM_ID",
        "OMIM_DiseaseName",
        "OMIM_Inheritance",
        "OMIM_Status",
    ]

    display_columns = [
        column
        for column in display_columns
        if column in dataframe.columns
    ]

    print("\nDatabase-integrated table:\n")
    print(
        dataframe[display_columns]
        .drop_duplicates()
        .to_string(index=False)
    )
    print(f"\nSaved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
