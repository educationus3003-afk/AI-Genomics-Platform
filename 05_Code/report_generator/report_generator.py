from __future__ import annotations

import html
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


BASE = Path.home() / "AI_Genomics_Platform"

INTERPRETATION_FILE = (
    BASE / "06_Results/Interpretation/ai_interpretations.json"
)

DATABASE_TABLE = (
    BASE / "06_Results/Master/database_integrated_table.csv"
)

OUTPUT_DIR = BASE / "08_Reports"
OUTPUT_HTML = OUTPUT_DIR / "genomic_ai_report.html"


def load_interpretations() -> list[dict[str, Any]]:
    if not INTERPRETATION_FILE.exists():
        raise FileNotFoundError(
            f"Interpretation file not found: {INTERPRETATION_FILE}"
        )

    with INTERPRETATION_FILE.open("r", encoding="utf-8") as file:
        return json.load(file)


def load_database_table() -> pd.DataFrame:
    if not DATABASE_TABLE.exists():
        raise FileNotFoundError(
            f"Database table not found: {DATABASE_TABLE}"
        )

    return pd.read_csv(DATABASE_TABLE)


def safe_text(value: Any) -> str:
    if value is None:
        return "Not available"

    if isinstance(value, float) and pd.isna(value):
        return "Not available"

    text = str(value).strip()
    return text if text else "Not available"


def evidence_row(
    label: str,
    value: Any,
) -> str:
    return (
        "<tr>"
        f"<th>{html.escape(label)}</th>"
        f"<td>{html.escape(safe_text(value))}</td>"
        "</tr>"
    )


def build_variant_section(
    interpretation: dict[str, Any],
    evidence: pd.Series | None,
) -> str:
    gene = safe_text(interpretation.get("gene"))
    disease = safe_text(interpretation.get("disease"))
    variant_id = safe_text(interpretation.get("variant_id"))
    classification = safe_text(
        interpretation.get("automated_classification")
    )

    if evidence is None:
        evidence = pd.Series(dtype=object)

    evidence_table = "".join(
        [
            evidence_row(
                "InterVar classification",
                evidence.get("INTERVAR_CLASSIFICATION"),
            ),
            evidence_row(
                "ACMG evidence",
                evidence.get("ACMG_EVIDENCE"),
            ),
            evidence_row(
                "ClinVar significance",
                evidence.get("ClinVar_ClinicalSignificance"),
            ),
            evidence_row(
                "ClinVar review status",
                evidence.get("ClinVar_ReviewStatus"),
            ),
            evidence_row(
                "ClinVar disease",
                evidence.get("ClinVar_Disease"),
            ),
            evidence_row(
                "gnomAD allele frequency",
                evidence.get("gnomAD_AlleleFrequency"),
            ),
            evidence_row(
                "AlphaMissense prediction",
                evidence.get("AlphaMissense_Prediction"),
            ),
            evidence_row(
                "REVEL score",
                evidence.get("REVEL_Score"),
            ),
            evidence_row(
                "CADD PHRED",
                evidence.get("CADD_PHRED"),
            ),
            evidence_row(
                "SpliceAI delta score",
                evidence.get("SpliceAI_DeltaScore"),
            ),
            evidence_row(
                "UniProt protein",
                evidence.get("UniProt_Protein"),
            ),
            evidence_row(
                "OMIM disease",
                evidence.get("OMIM_DiseaseName"),
            ),
        ]
    )

    clinical_features = interpretation.get(
        "clinical_interpretation",
        "Not available",
    )

    return f"""
    <section class="variant-card">
        <div class="variant-header">
            <div>
                <h2>{html.escape(gene)} — {html.escape(disease)}</h2>
                <p class="variant-id">{html.escape(variant_id)}</p>
            </div>
            <span class="classification">
                {html.escape(classification)}
            </span>
        </div>

        <div class="summary-grid">
            <div class="summary-box">
                <h3>Gene summary</h3>
                <p>{html.escape(
                    safe_text(interpretation.get("gene_summary"))
                )}</p>
            </div>

            <div class="summary-box">
                <h3>Variant summary</h3>
                <p>{html.escape(
                    safe_text(interpretation.get("variant_summary"))
                )}</p>
            </div>

            <div class="summary-box">
                <h3>Disease summary</h3>
                <p>{html.escape(
                    safe_text(interpretation.get("disease_summary"))
                )}</p>
            </div>

            <div class="summary-box">
                <h3>Patient-friendly explanation</h3>
                <p>{html.escape(
                    safe_text(
                        interpretation.get(
                            "patient_friendly_explanation"
                        )
                    )
                )}</p>
            </div>
        </div>

        <h3>Integrated evidence</h3>
        <table class="evidence-table">
            <tbody>
                {evidence_table}
            </tbody>
        </table>

        <h3>Clinical interpretation</h3>
        <p>{html.escape(safe_text(clinical_features))}</p>

        <div class="warning">
            This automated interpretation is for research and
            decision-support purposes only. It does not replace
            review by a qualified clinical genetics professional.
        </div>
    </section>
    """


def build_report() -> str:
    interpretations = load_interpretations()
    database_table = load_database_table()

    variant_sections: list[str] = []

    for interpretation in interpretations:
        variant_id = safe_text(interpretation.get("variant_id"))

        matches = database_table[
            database_table["VARIANT_ID"].astype(str) == variant_id
        ]

        evidence = None
        if not matches.empty:
            evidence = matches.iloc[0]

        variant_sections.append(
            build_variant_section(
                interpretation=interpretation,
                evidence=evidence,
            )
        )

    generated_time = datetime.now().strftime(
        "%Y-%m-%d %H:%M:%S"
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport"
          content="width=device-width, initial-scale=1.0">
    <title>AI Genomics Interpretation Report</title>

    <style>
        body {{
            font-family: Arial, Helvetica, sans-serif;
            margin: 0;
            background: #f4f7fb;
            color: #1f2937;
            line-height: 1.6;
        }}

        .container {{
            max-width: 1100px;
            margin: 0 auto;
            padding: 32px;
        }}

        .report-header {{
            background: #ffffff;
            border-radius: 14px;
            padding: 30px;
            margin-bottom: 24px;
            box-shadow: 0 4px 18px rgba(0, 0, 0, 0.08);
        }}

        .report-header h1 {{
            margin: 0 0 10px;
            font-size: 32px;
        }}

        .subtitle {{
            color: #4b5563;
            margin-bottom: 14px;
        }}

        .metadata {{
            font-size: 14px;
            color: #6b7280;
        }}

        .variant-card {{
            background: #ffffff;
            border-radius: 14px;
            padding: 26px;
            margin-bottom: 26px;
            box-shadow: 0 4px 18px rgba(0, 0, 0, 0.08);
        }}

        .variant-header {{
            display: flex;
            justify-content: space-between;
            gap: 20px;
            align-items: flex-start;
            border-bottom: 1px solid #e5e7eb;
            padding-bottom: 16px;
            margin-bottom: 20px;
        }}

        .variant-header h2 {{
            margin: 0;
            font-size: 24px;
        }}

        .variant-id {{
            margin: 5px 0 0;
            font-family: monospace;
            color: #4b5563;
        }}

        .classification {{
            display: inline-block;
            padding: 8px 12px;
            border-radius: 999px;
            background: #e8eefc;
            font-weight: bold;
            white-space: nowrap;
        }}

        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 16px;
            margin-bottom: 22px;
        }}

        .summary-box {{
            background: #f8fafc;
            border: 1px solid #e5e7eb;
            border-radius: 10px;
            padding: 16px;
        }}

        .summary-box h3 {{
            margin-top: 0;
        }}

        .evidence-table {{
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 20px;
        }}

        .evidence-table th,
        .evidence-table td {{
            border: 1px solid #d1d5db;
            text-align: left;
            vertical-align: top;
            padding: 10px;
        }}

        .evidence-table th {{
            width: 34%;
            background: #f3f4f6;
        }}

        .warning {{
            margin-top: 18px;
            background: #fff7ed;
            border-left: 5px solid #f97316;
            padding: 14px;
            border-radius: 6px;
            font-size: 14px;
        }}

        .footer {{
            text-align: center;
            color: #6b7280;
            font-size: 13px;
            margin-top: 25px;
        }}

        @media (max-width: 760px) {{
            .summary-grid {{
                grid-template-columns: 1fr;
            }}

            .variant-header {{
                flex-direction: column;
            }}
        }}

        @media print {{
            body {{
                background: #ffffff;
            }}

            .variant-card,
            .report-header {{
                box-shadow: none;
                border: 1px solid #d1d5db;
            }}

            .variant-card {{
                page-break-inside: avoid;
            }}
        }}
    </style>
</head>

<body>
    <main class="container">
        <header class="report-header">
            <h1>AI Genomics Interpretation Report</h1>
            <p class="subtitle">
                Automated evidence integration, disease-context
                interpretation, and clinician-oriented reporting
                for genomic variants.
            </p>
            <p class="metadata">
                Generated: {html.escape(generated_time)}<br>
                Variants analysed: {len(interpretations)}
            </p>
        </header>

        {''.join(variant_sections)}

        <footer class="footer">
            AI Genomics Platform — Research Prototype
        </footer>
    </main>
</body>
</html>
"""


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    report = build_report()

    OUTPUT_HTML.write_text(
        report,
        encoding="utf-8",
    )

    print(f"\nReport generated successfully:")
    print(OUTPUT_HTML)


if __name__ == "__main__":
    main()
