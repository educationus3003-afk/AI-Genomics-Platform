from pathlib import Path
import pandas as pd


BASE = Path.home() / "AI_Genomics_Platform"

MASTER_INPUT = BASE / "06_Results/Master/master_annotation_table.csv"
OUTPUT_FILE = BASE / "06_Results/Master/standardized_annotation_table.csv"


def load_master_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    df = pd.read_csv(path)

    required_columns = ["CHROM", "POS", "REF", "ALT", "GENE"]
    missing = [column for column in required_columns if column not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    return df


def standardize_annotations(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    df["CHROM"] = df["CHROM"].astype(str)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["REF"] = df["REF"].astype(str)
    df["ALT"] = df["ALT"].astype(str)
    df["GENE"] = df["GENE"].astype(str).str.upper()

    df["VARIANT_ID"] = (
        "chr"
        + df["CHROM"].str.replace("chr", "", regex=False)
        + ":"
        + df["POS"].astype("Int64").astype(str)
        + ":"
        + df["REF"]
        + ">"
        + df["ALT"]
    )

    if "INTERVAR_CLASSIFICATION" not in df.columns:
        df["INTERVAR_CLASSIFICATION"] = "Not available"

    if "ACMG_EVIDENCE" not in df.columns:
        df["ACMG_EVIDENCE"] = "Not available"

    if "SOURCE" not in df.columns:
        df["SOURCE"] = "Not available"

    df["INTERVAR_CLASSIFICATION"] = (
        df["INTERVAR_CLASSIFICATION"]
        .fillna("Not available")
        .replace("Pending", "Pending web review")
    )

    df["ACMG_EVIDENCE"] = (
        df["ACMG_EVIDENCE"]
        .fillna("Not available")
        .replace("Pending", "Pending web review")
    )

    ordered_columns = [
        "VARIANT_ID",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "GENE",
        "INTERVAR_CLASSIFICATION",
        "ACMG_EVIDENCE",
        "SOURCE",
    ]

    return df[ordered_columns]


def main() -> None:
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    master = load_master_table(MASTER_INPUT)
    standardized = standardize_annotations(master)
    standardized.to_csv(OUTPUT_FILE, index=False)

    print("\nStandardized annotation table:")
    print(standardized.to_string(index=False))
    print(f"\nSaved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
