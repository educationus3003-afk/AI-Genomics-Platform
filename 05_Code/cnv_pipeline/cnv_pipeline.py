from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd


BASE = Path.home() / "AI_Genomics_Platform"
DEFAULT_INPUT = BASE / "04_Data/Test/CNV/sample_cnv.bed"
DEFAULT_OUTPUT = BASE / "06_Results/CNV"

REQUIRED_COLUMNS = [
    "CHROM",
    "START",
    "END",
    "CNV_TYPE",
    "SAMPLE",
]


def executable_status(names: list[str]) -> tuple[bool, str | None]:
    for name in names:
        path = shutil.which(name)
        if path:
            return True, path
    return False, None


def validate_cnv_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"CNV input file not found: {path}")

    dataframe = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        dtype={
            "CHROM": str,
            "CNV_TYPE": str,
            "SAMPLE": str,
        },
    )

    missing = [
        column
        for column in REQUIRED_COLUMNS
        if column not in dataframe.columns
    ]

    if missing:
        raise ValueError(
            f"CNV input is missing required columns: {missing}"
        )

    dataframe["START"] = pd.to_numeric(
        dataframe["START"],
        errors="raise",
    ).astype(int)

    dataframe["END"] = pd.to_numeric(
        dataframe["END"],
        errors="raise",
    ).astype(int)

    invalid_coordinates = dataframe["END"] <= dataframe["START"]

    if invalid_coordinates.any():
        rows = list(dataframe.index[invalid_coordinates] + 2)
        raise ValueError(
            f"END must be greater than START. Invalid rows: {rows}"
        )

    dataframe["CNV_TYPE"] = (
        dataframe["CNV_TYPE"]
        .str.upper()
        .str.strip()
    )

    accepted = {"DEL", "DUP"}

    invalid_types = ~dataframe["CNV_TYPE"].isin(accepted)

    if invalid_types.any():
        bad = sorted(
            dataframe.loc[
                invalid_types,
                "CNV_TYPE",
            ].unique()
        )
        raise ValueError(
            f"Unsupported CNV types: {bad}. Use DEL or DUP."
        )

    dataframe["CHROM"] = (
        dataframe["CHROM"]
        .astype(str)
        .str.replace("chr", "", regex=False)
    )

    dataframe["CNV_SIZE_BP"] = (
        dataframe["END"] - dataframe["START"]
    )

    dataframe["CNV_ID"] = (
        "chr"
        + dataframe["CHROM"]
        + ":"
        + dataframe["START"].astype(str)
        + "-"
        + dataframe["END"].astype(str)
        + ":"
        + dataframe["CNV_TYPE"]
    )

    return dataframe


def get_tool_status() -> list[dict[str, Any]]:
    definitions = [
        {
            "tool": "AnnotSV",
            "executables": ["AnnotSV", "annotSV"],
            "purpose": "Structural-variant annotation",
            "expected_input": "BED or SV-VCF",
            "expected_output": "Annotated TSV",
        },
        {
            "tool": "ClassifyCNV",
            "executables": ["ClassifyCNV", "classifycnv"],
            "purpose": "ACMG/ClinGen CNV classification",
            "expected_input": "BED",
            "expected_output": "Classification and evidence table",
        },
        {
            "tool": "ISV-CNV",
            "executables": ["isvcnv", "ISV-CNV"],
            "purpose": "CNV interpretation support",
            "expected_input": "CNV coordinates and annotations",
            "expected_output": "Interpretation evidence",
        },
        {
            "tool": "ClinGen CNV",
            "executables": [],
            "purpose": "Dosage sensitivity and CNV evidence",
            "expected_input": "CNV region",
            "expected_output": "Manual or API-derived evidence",
        },
    ]

    statuses: list[dict[str, Any]] = []

    for definition in definitions:
        executables = definition["executables"]

        if executables:
            installed, path = executable_status(executables)
        else:
            installed, path = False, None

        status = {
            **definition,
            "installed": installed,
            "executable_path": path or "Not detected",
            "integration_status": (
                "Executable detected"
                if installed
                else "Documented integration; external setup required"
            ),
        }

        statuses.append(status)

    return statuses


def add_interpretation_columns(
    dataframe: pd.DataFrame,
) -> pd.DataFrame:
    result = dataframe.copy()

    result["AnnotSV_Status"] = "Not executed"
    result["ClassifyCNV_Status"] = "Not executed"
    result["ISV_CNV_Status"] = "Not executed"
    result["ClinGen_Status"] = "Manual review required"

    result["Automated_CNV_Classification"] = (
        "Unclassified — external evidence required"
    )

    result["Interpretation_Note"] = (
        "CNV requires gene-content, dosage-sensitivity, "
        "population-frequency and phenotype review."
    )

    return result


def write_html_report(
    cnvs: pd.DataFrame,
    tool_status: pd.DataFrame,
    output_file: Path,
) -> None:
    cnv_table = cnvs.to_html(
        index=False,
        classes="data-table",
        border=0,
    )

    status_table = tool_status.to_html(
        index=False,
        classes="data-table",
        border=0,
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>GeneInsightAI CNV Report</title>
<style>
body {{
    font-family: Arial, sans-serif;
    margin: 40px;
    background: #f5f7fa;
    color: #1f2937;
}}
.container {{
    max-width: 1500px;
    margin: auto;
    background: white;
    padding: 30px;
    border-radius: 12px;
}}
h1, h2 {{
    color: #17365d;
}}
.notice {{
    padding: 14px;
    background: #fff4ce;
    border-left: 5px solid #d69e00;
    margin-bottom: 24px;
}}
.data-table {{
    border-collapse: collapse;
    width: 100%;
    margin-bottom: 30px;
    font-size: 13px;
}}
.data-table th,
.data-table td {{
    border: 1px solid #d1d5db;
    padding: 8px;
    text-align: left;
}}
.data-table th {{
    background: #17365d;
    color: white;
}}
.data-table tr:nth-child(even) {{
    background: #f3f4f6;
}}
.footer {{
    margin-top: 30px;
    font-size: 12px;
    color: #6b7280;
}}
</style>
</head>
<body>
<div class="container">

<h1>GeneInsightAI CNV Interpretation Report</h1>

<div class="notice">
<strong>Research-use notice:</strong>
The workflow records tool availability and prepares CNVs for
interpretation. A qualified genetics professional must review
all classifications.
</div>

<h2>CNV Input Summary</h2>
<p>Total CNVs: <strong>{len(cnvs)}</strong></p>
{cnv_table}

<h2>CNV Tool Integration Status</h2>
{status_table}

<h2>Workflow</h2>
<p>
CNV BED/SV-VCF → AnnotSV → ClassifyCNV → ISV-CNV →
ClinGen evidence → CNV interpretation report
</p>

<div class="footer">
Generated by GeneInsightAI CNV branch.
Unavailable external tools are reported honestly and no
pathogenicity evidence is fabricated.
</div>

</div>
</body>
</html>
"""

    output_file.write_text(
        html,
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="GeneInsightAI CNV analysis framework"
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Input tab-separated CNV BED-style file",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output directory",
    )

    args = parser.parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("GENEINSIGHTAI CNV PIPELINE")
    print("=" * 70)

    cnvs = validate_cnv_table(args.input)
    cnvs = add_interpretation_columns(cnvs)

    statuses = get_tool_status()
    tool_status = pd.DataFrame(statuses)

    cnv_csv = output_dir / "cnv_interpretation.csv"
    tool_csv = output_dir / "cnv_tool_status.csv"
    json_file = output_dir / "cnv_interpretation.json"
    html_file = output_dir / "cnv_report.html"

    cnvs.to_csv(cnv_csv, index=False)
    tool_status.to_csv(tool_csv, index=False)

    payload = {
        "input_file": str(args.input),
        "cnv_count": int(len(cnvs)),
        "cnvs": cnvs.to_dict(orient="records"),
        "tools": statuses,
        "warning": (
            "Research-use output. CNV interpretation requires "
            "professional review."
        ),
    }

    json_file.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )

    write_html_report(
        cnvs=cnvs,
        tool_status=tool_status,
        output_file=html_file,
    )

    print(f"\nValidated CNVs: {len(cnvs)}")
    print("\nTool status:")

    for record in statuses:
        symbol = "READY" if record["installed"] else "SETUP REQUIRED"
        print(f"  [{symbol}] {record['tool']}")

    print("\nGenerated outputs:")
    print(f"  [READY] {cnv_csv}")
    print(f"  [READY] {tool_csv}")
    print(f"  [READY] {json_file}")
    print(f"  [READY] {html_file}")

    print("\nCNV pipeline completed successfully.")


if __name__ == "__main__":
    main()
