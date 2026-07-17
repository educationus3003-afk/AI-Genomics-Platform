from __future__ import annotations

import base64
import html
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd


BASE = Path.home() / "AI_Genomics_Platform"

OUTPUT_FILE = (
    BASE
    / "08_Reports"
    / "AI_Genomics_Platform_Complete_Report.html"
)

SUBMISSION_DIR = (
    BASE
    / "09_Presentation"
    / "Final_Submission"
)

VERIFIED_TABLE = (
    SUBMISSION_DIR
    / "final_verified_variant_table.csv"
)

MASTER_TABLE = (
    BASE
    / "06_Results"
    / "Master"
    / "database_integrated_table.csv"
)

AI_INTERPRETATIONS = (
    BASE
    / "06_Results"
    / "Interpretation"
    / "ai_interpretations.csv"
)

CNV_TABLE = (
    BASE
    / "06_Results"
    / "CNV"
    / "cnv_interpretation.csv"
)

CNV_STATUS = (
    BASE
    / "06_Results"
    / "CNV"
    / "cnv_tool_status.csv"
)

COMPARISON_DIR = (
    BASE
    / "06_Results"
    / "Comparison"
)

FIGURE_FILES = {
    "Gene Distribution": (
        BASE
        / "07_Visualizations"
        / "gene_distribution.png"
    ),
    "Classification Distribution": (
        BASE
        / "07_Visualizations"
        / "classification_pie.png"
    ),
    "Database Coverage": (
        BASE
        / "07_Visualizations"
        / "database_coverage.png"
    ),
}

DISEASE_DETAILS = {
    "CFTR": {
        "disease": "Cystic fibrosis",
        "inheritance": "Autosomal recessive",
        "overview": (
            "Cystic fibrosis is caused by pathogenic variants in "
            "CFTR, which encodes an epithelial chloride and "
            "bicarbonate channel. Dysfunction leads to thick "
            "secretions affecting the lungs, pancreas and other "
            "organs."
        ),
        "mechanism": (
            "Loss of CFTR channel activity disrupts epithelial "
            "ion transport and fluid balance."
        ),
    },
    "HBB": {
        "disease": "Sickle cell disease",
        "inheritance": "Autosomal recessive",
        "overview": (
            "Sickle cell disease is caused by pathogenic HBB "
            "variants that alter beta-globin structure and "
            "haemoglobin behaviour."
        ),
        "mechanism": (
            "Abnormal haemoglobin polymerisation causes red-cell "
            "deformation, haemolysis and vaso-occlusion."
        ),
    },
    "MECP2": {
        "disease": "Rett syndrome",
        "inheritance": "X-linked dominant",
        "overview": (
            "Rett syndrome is a neurodevelopmental disorder "
            "associated mainly with pathogenic MECP2 variants."
        ),
        "mechanism": (
            "MECP2 dysfunction disturbs transcriptional regulation "
            "and neuronal maturation."
        ),
    },
    "FBN1": {
        "disease": "Marfan syndrome",
        "inheritance": "Autosomal dominant",
        "overview": (
            "Marfan syndrome is a connective-tissue disorder "
            "caused by pathogenic variants in FBN1."
        ),
        "mechanism": (
            "Fibrillin-1 dysfunction affects extracellular "
            "microfibrils and TGF-beta signalling."
        ),
    },
}


def safe_read_csv(
    path: Path,
    *,
    sep: str = ",",
) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()

    try:
        return pd.read_csv(path, sep=sep)
    except Exception as error:
        print(f"Warning: could not read {path}: {error}")
        return pd.DataFrame()


def normalise_missing(value: object) -> str:
    if pd.isna(value):
        return "Unavailable"

    text = str(value).strip()

    if text.lower() in {
        "",
        "nan",
        "none",
        "not retrieved",
        "not available",
    }:
        return "Unavailable"

    return text


def first_present(
    row: pd.Series,
    columns: Iterable[str],
    default: str = "Unavailable",
) -> str:
    for column in columns:
        if column in row.index:
            value = normalise_missing(row[column])
            if value != "Unavailable":
                return value

    return default


def image_data_uri(path: Path) -> str | None:
    if not path.exists():
        return None

    encoded = base64.b64encode(
        path.read_bytes()
    ).decode("ascii")

    suffix = path.suffix.lower().replace(".", "")
    mime = "jpeg" if suffix in {"jpg", "jpeg"} else suffix

    return f"data:image/{mime};base64,{encoded}"


def dataframe_html(
    dataframe: pd.DataFrame,
    *,
    empty_message: str,
    max_rows: int | None = None,
) -> str:
    if dataframe.empty:
        return (
            '<div class="empty-note">'
            f"{html.escape(empty_message)}"
            "</div>"
        )

    view = dataframe.copy()

    if max_rows is not None:
        view = view.head(max_rows)

    view = view.map(normalise_missing)

    return view.to_html(
        index=False,
        classes="report-table",
        border=0,
        escape=True,
    )


def load_variant_table() -> pd.DataFrame:
    if VERIFIED_TABLE.exists():
        dataframe = safe_read_csv(VERIFIED_TABLE)
    else:
        dataframe = safe_read_csv(MASTER_TABLE)

    if dataframe.empty:
        return dataframe

    if "GENE" in dataframe.columns:
        dataframe["GENE"] = (
            dataframe["GENE"]
            .astype(str)
            .str.upper()
            .str.strip()
        )

    return dataframe


def make_metric_cards(
    variants: pd.DataFrame,
) -> str:
    gene_count = (
        variants["GENE"].nunique()
        if "GENE" in variants.columns
        else 0
    )

    variant_count = len(variants)

    available_omim = 0
    if "OMIM_OMIM_ID" in variants.columns:
        available_omim = sum(
            normalise_missing(value) != "Unavailable"
            for value in variants["OMIM_OMIM_ID"]
        )

    return f"""
    <div class="metrics">
        <div class="metric">
            <div class="metric-value">{variant_count}</div>
            <div class="metric-label">Validated variants</div>
        </div>

        <div class="metric">
            <div class="metric-value">{gene_count}</div>
            <div class="metric-label">Disease genes</div>
        </div>

        <div class="metric">
            <div class="metric-value">{available_omim}</div>
            <div class="metric-label">OMIM associations</div>
        </div>

        <div class="metric">
            <div class="metric-value">2</div>
            <div class="metric-label">Analysis branches</div>
        </div>
    </div>
    """


def make_disease_sections(
    variants: pd.DataFrame,
    interpretations: pd.DataFrame,
) -> str:
    sections: list[str] = []

    for gene, details in DISEASE_DETAILS.items():
        gene_rows = pd.DataFrame()

        if (
            not variants.empty
            and "GENE" in variants.columns
        ):
            gene_rows = variants[
                variants["GENE"] == gene
            ]

        if gene_rows.empty:
            row = pd.Series(dtype=object)
        else:
            row = gene_rows.iloc[0]

        variant_id = first_present(
            row,
            ["VARIANT_ID"],
        )

        omim_id = first_present(
            row,
            ["OMIM_OMIM_ID", "OMIM_ID"],
        )

        disease_name = first_present(
            row,
            ["OMIM_DiseaseName", "Disease"],
            details["disease"],
        )

        inheritance = first_present(
            row,
            ["OMIM_Inheritance", "Inheritance"],
            details["inheritance"],
        )

        consequence = first_present(
            row,
            [
                "Final_Consequence",
                "VEP_Consequence",
                "Input_Consequence",
            ],
        )

        hgvsc = first_present(
            row,
            [
                "Final_HGVSc",
                "VEP_HGVSc",
                "Input_HGVSc",
            ],
        )

        hgvsp = first_present(
            row,
            [
                "Final_HGVSp",
                "VEP_HGVSp",
                "Input_HGVSp",
            ],
        )

        dbsnp = first_present(
            row,
            [
                "Final_dbSNP_ID",
                "dbSNP_ID",
                "Input_dbSNP_ID",
            ],
        )

        clinvar = first_present(
            row,
            [
                "ClinVar_ClinicalSignificance",
                "Input_ClinVar",
            ],
        )

        gnomad = first_present(
            row,
            [
                "gnomAD_AlleleFrequency",
                "gnomAD_AF",
                "Input_gnomAD_AF",
            ],
        )

        alpha = first_present(
            row,
            [
                "AlphaMissense_Prediction",
                "AlphaMissense_Score",
                "Input_AlphaMissense",
            ],
        )

        annotation_source = first_present(
            row,
            [
                "Annotation_Source",
                "VEP_Source",
                "SOURCE",
            ],
            (
                "Validated input metadata and "
                "available annotation outputs"
            ),
        )

        ai_text = (
            "Interpretation was generated through the "
            "platform's AI-assisted reporting module."
        )

        if not interpretations.empty:
            matching = pd.DataFrame()

            if "GENE" in interpretations.columns:
                matching = interpretations[
                    interpretations["GENE"]
                    .astype(str)
                    .str.upper()
                    == gene
                ]

            if not matching.empty:
                ai_candidates = [
                    "AI_Interpretation",
                    "Interpretation",
                    "Clinical_Interpretation",
                    "Summary",
                    "AI_Summary",
                ]

                ai_text = first_present(
                    matching.iloc[0],
                    ai_candidates,
                    ai_text,
                )

        evidence_rows = [
            ("Variant", variant_id),
            ("Consequence", consequence),
            ("HGVSc", hgvsc),
            ("HGVSp", hgvsp),
            ("dbSNP", dbsnp),
            ("ClinVar", clinvar),
            ("gnomAD allele frequency", gnomad),
            ("AlphaMissense", alpha),
            ("OMIM ID", omim_id),
            ("Inheritance", inheritance),
            ("Annotation source", annotation_source),
        ]

        evidence_html = "".join(
            f"""
            <tr>
                <th>{html.escape(label)}</th>
                <td>{html.escape(normalise_missing(value))}</td>
            </tr>
            """
            for label, value in evidence_rows
        )

        sections.append(
            f"""
            <section class="disease-section">
                <div class="disease-heading">
                    <div>
                        <div class="gene-badge">{gene}</div>
                        <h3>{html.escape(disease_name)}</h3>
                    </div>
                    <span class="inheritance-tag">
                        {html.escape(inheritance)}
                    </span>
                </div>

                <div class="summary-grid">
                    <article class="summary-card">
                        <h4>Disease overview</h4>
                        <p>{html.escape(details["overview"])}</p>
                    </article>

                    <article class="summary-card">
                        <h4>Molecular mechanism</h4>
                        <p>{html.escape(details["mechanism"])}</p>
                    </article>
                </div>

                <h4>Integrated evidence</h4>

                <table class="evidence-table">
                    {evidence_html}
                </table>

                <div class="interpretation-box">
                    <h4>AI-assisted interpretation</h4>
                    <p>{html.escape(ai_text)}</p>
                </div>
            </section>
            """
        )

    return "\n".join(sections)


def make_figures_section() -> str:
    figures: list[str] = []

    for title, path in FIGURE_FILES.items():
        uri = image_data_uri(path)

        if uri is None:
            continue

        figures.append(
            f"""
            <figure class="figure-card">
                <img
                    src="{uri}"
                    alt="{html.escape(title)}"
                >
                <figcaption>
                    {html.escape(title)}
                </figcaption>
            </figure>
            """
        )

    if not figures:
        return (
            '<div class="empty-note">'
            "Visualization figures were not available."
            "</div>"
        )

    return (
        '<div class="figure-grid">'
        + "".join(figures)
        + "</div>"
    )


def make_comparison_section() -> str:
    tables: list[str] = []

    comparison_files = [
        (
            "Normalization comparison",
            COMPARISON_DIR / "normalization_comparison.tsv",
        ),
        (
            "Reference validation",
            COMPARISON_DIR / "reference_check.tsv",
        ),
        (
            "VEP annotations",
            COMPARISON_DIR / "vep_annotations.tsv",
        ),
        (
            "SnpEff annotations",
            COMPARISON_DIR / "snpeff_annotations.tsv",
        ),
    ]

    for title, path in comparison_files:
        dataframe = safe_read_csv(path, sep="\t")

        if dataframe.empty:
            continue

        tables.append(
            f"""
            <div class="comparison-block">
                <h3>{html.escape(title)}</h3>
                {dataframe_html(
                    dataframe,
                    empty_message="No comparison data available.",
                    max_rows=20,
                )}
            </div>
            """
        )

    if not tables:
        return (
            '<div class="empty-note">'
            "No comparison outputs were available."
            "</div>"
        )

    return "\n".join(tables)


def make_tools_table() -> str:
    tools = pd.DataFrame(
        [
            {
                "Tool / Database": "bcftools",
                "Contribution": (
                    "VCF normalization and reference checking"
                ),
                "Status": "Implemented",
            },
            {
                "Tool / Database": "VEP",
                "Contribution": (
                    "Transcript consequence and HGVS annotation"
                ),
                "Status": (
                    "Integration implemented; local cache "
                    "required for a fresh offline run"
                ),
            },
            {
                "Tool / Database": "SnpEff",
                "Contribution": (
                    "Functional consequence comparison"
                ),
                "Status": "Implemented",
            },
            {
                "Tool / Database": "ClinVar",
                "Contribution": (
                    "Clinical significance and review evidence"
                ),
                "Status": (
                    "Framework implemented; availability varies"
                ),
            },
            {
                "Tool / Database": "gnomAD",
                "Contribution": (
                    "Population allele frequency evidence"
                ),
                "Status": (
                    "Framework implemented; availability varies"
                ),
            },
            {
                "Tool / Database": "dbSNP",
                "Contribution": (
                    "Known-variant identifier"
                ),
                "Status": "Integrated where available",
            },
            {
                "Tool / Database": "OMIM",
                "Contribution": (
                    "Gene-disease and inheritance reference data"
                ),
                "Status": "Local reference metadata integrated",
            },
            {
                "Tool / Database": "UniProt",
                "Contribution": (
                    "Protein name and functional context"
                ),
                "Status": "Integration framework implemented",
            },
            {
                "Tool / Database": "AlphaMissense",
                "Contribution": (
                    "Missense pathogenicity prediction"
                ),
                "Status": (
                    "Applicable only to eligible missense variants"
                ),
            },
            {
                "Tool / Database": "AnnotSV",
                "Contribution": (
                    "Structural-variant annotation"
                ),
                "Status": "External setup required",
            },
            {
                "Tool / Database": "ClassifyCNV",
                "Contribution": (
                    "ACMG/ClinGen CNV classification"
                ),
                "Status": "External setup required",
            },
            {
                "Tool / Database": "ISV-CNV",
                "Contribution": "CNV interpretation support",
                "Status": "External setup required",
            },
            {
                "Tool / Database": "ClinGen CNV",
                "Contribution": (
                    "Dosage-sensitivity and CNV evidence"
                ),
                "Status": "Manual/external evidence required",
            },
        ]
    )

    return dataframe_html(
        tools,
        empty_message="No tool information available.",
    )


def build_report() -> str:
    variants = load_variant_table()
    interpretations = safe_read_csv(AI_INTERPRETATIONS)
    cnvs = safe_read_csv(CNV_TABLE)
    cnv_status = safe_read_csv(CNV_STATUS)

    variant_columns = [
        "VARIANT_ID",
        "GENE",
        "OMIM_OMIM_ID",
        "OMIM_DiseaseName",
        "OMIM_Inheritance",
        "Final_dbSNP_ID",
        "Final_Consequence",
        "Final_HGVSc",
        "Final_HGVSp",
        "Annotation_Source",
    ]

    available_variant_columns = [
        column
        for column in variant_columns
        if column in variants.columns
    ]

    variant_summary = (
        variants[available_variant_columns]
        if available_variant_columns
        else variants
    )

    generated = datetime.now().strftime(
        "%d %B %Y, %H:%M"
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta
    name="viewport"
    content="width=device-width, initial-scale=1.0"
>
<title>
AI Genomics Platform — Complete Interpretation Report
</title>

<style>
:root {{
    --navy: #102a43;
    --blue: #1f5f99;
    --light-blue: #eaf4fb;
    --green: #1f7a55;
    --orange: #b45309;
    --red: #b42318;
    --grey: #5f6b76;
    --light-grey: #f4f7fa;
    --border: #d8e0e8;
    --white: #ffffff;
}}

* {{
    box-sizing: border-box;
}}

body {{
    margin: 0;
    background: #eef2f6;
    color: #1f2933;
    font-family:
        Inter,
        "Segoe UI",
        Arial,
        sans-serif;
    line-height: 1.6;
}}

.report-container {{
    max-width: 1450px;
    margin: 28px auto;
    background: var(--white);
    box-shadow: 0 10px 35px rgba(16, 42, 67, 0.12);
}}

.cover {{
    background:
        linear-gradient(
            135deg,
            #102a43 0%,
            #1f5f99 68%,
            #2f80b9 100%
        );
    color: white;
    padding: 70px 70px 58px;
}}

.cover-label {{
    display: inline-block;
    margin-bottom: 22px;
    padding: 7px 14px;
    border: 1px solid rgba(255,255,255,0.4);
    border-radius: 20px;
    font-size: 13px;
    font-weight: 700;
    letter-spacing: 1.2px;
    text-transform: uppercase;
}}

.cover h1 {{
    max-width: 950px;
    margin: 0 0 18px;
    font-size: 46px;
    line-height: 1.12;
}}

.cover-subtitle {{
    max-width: 950px;
    margin: 0;
    font-size: 20px;
    color: #dbeaf5;
}}

.cover-meta {{
    margin-top: 36px;
    color: #dbeaf5;
    font-size: 14px;
}}

.content {{
    padding: 45px 65px 70px;
}}

section {{
    margin-bottom: 55px;
}}

.section-title {{
    margin: 0 0 22px;
    padding-bottom: 10px;
    border-bottom: 3px solid var(--blue);
    color: var(--navy);
    font-size: 29px;
}}

.lead {{
    max-width: 1100px;
    font-size: 17px;
    color: #425466;
}}

.metrics {{
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 18px;
    margin-top: 28px;
}}

.metric {{
    padding: 25px;
    border: 1px solid var(--border);
    border-radius: 12px;
    background: var(--light-grey);
}}

.metric-value {{
    color: var(--blue);
    font-size: 34px;
    font-weight: 800;
}}

.metric-label {{
    color: var(--grey);
    font-size: 14px;
    font-weight: 600;
}}

.workflow {{
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 10px;
    margin-top: 26px;
}}

.workflow-step {{
    padding: 13px 18px;
    border-radius: 8px;
    background: var(--light-blue);
    color: var(--navy);
    font-weight: 700;
}}

.workflow-arrow {{
    color: var(--blue);
    font-size: 24px;
    font-weight: 800;
}}

.report-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 13px;
    overflow-wrap: anywhere;
}}

.report-table th {{
    padding: 11px 10px;
    background: var(--navy);
    color: white;
    text-align: left;
    vertical-align: top;
}}

.report-table td {{
    padding: 10px;
    border-bottom: 1px solid var(--border);
    vertical-align: top;
}}

.report-table tr:nth-child(even) {{
    background: #f7f9fb;
}}

.disease-section {{
    padding: 26px;
    border: 1px solid var(--border);
    border-radius: 14px;
    background: #ffffff;
    box-shadow: 0 5px 18px rgba(16, 42, 67, 0.06);
}}

.disease-heading {{
    display: flex;
    justify-content: space-between;
    align-items: start;
    gap: 18px;
    margin-bottom: 22px;
}}

.disease-heading h3 {{
    margin: 8px 0 0;
    color: var(--navy);
    font-size: 25px;
}}

.gene-badge {{
    display: inline-block;
    padding: 5px 11px;
    border-radius: 6px;
    background: var(--blue);
    color: white;
    font-weight: 800;
    letter-spacing: 0.7px;
}}

.inheritance-tag {{
    padding: 7px 12px;
    border-radius: 18px;
    background: #ecfdf3;
    color: var(--green);
    font-size: 13px;
    font-weight: 700;
}}

.summary-grid {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 17px;
    margin-bottom: 24px;
}}

.summary-card {{
    padding: 17px;
    border-left: 4px solid var(--blue);
    background: var(--light-grey);
}}

.summary-card h4 {{
    margin: 0 0 7px;
    color: var(--navy);
}}

.summary-card p {{
    margin: 0;
}}

.evidence-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 14px;
}}

.evidence-table th {{
    width: 240px;
    padding: 10px;
    background: #f1f5f9;
    color: var(--navy);
    text-align: left;
    border: 1px solid var(--border);
}}

.evidence-table td {{
    padding: 10px;
    border: 1px solid var(--border);
    overflow-wrap: anywhere;
}}

.interpretation-box {{
    margin-top: 20px;
    padding: 18px 20px;
    border-left: 5px solid var(--green);
    background: #f0fdf7;
}}

.interpretation-box h4 {{
    margin: 0 0 7px;
    color: var(--green);
}}

.interpretation-box p {{
    margin: 0;
}}

.figure-grid {{
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 20px;
}}

.figure-card {{
    margin: 0;
    padding: 14px;
    border: 1px solid var(--border);
    border-radius: 12px;
    background: white;
}}

.figure-card img {{
    display: block;
    width: 100%;
    height: auto;
}}

.figure-card figcaption {{
    margin-top: 10px;
    color: var(--navy);
    text-align: center;
    font-weight: 700;
}}

.comparison-block {{
    margin-bottom: 30px;
}}

.comparison-block h3 {{
    color: var(--navy);
}}

.status-note {{
    padding: 18px;
    border-left: 5px solid var(--orange);
    background: #fff7ed;
}}

.status-note strong {{
    color: var(--orange);
}}

.warning {{
    padding: 20px;
    border-left: 5px solid var(--red);
    background: #fff1f0;
}}

.empty-note {{
    padding: 15px;
    border: 1px dashed var(--border);
    background: var(--light-grey);
    color: var(--grey);
}}

.footer {{
    padding: 24px 65px;
    background: var(--navy);
    color: #dbeaf5;
    font-size: 13px;
}}

@media (max-width: 900px) {{
    .cover {{
        padding: 45px 28px;
    }}

    .cover h1 {{
        font-size: 35px;
    }}

    .content {{
        padding: 35px 24px;
    }}

    .metrics,
    .summary-grid,
    .figure-grid {{
        grid-template-columns: 1fr;
    }}

    .disease-heading {{
        flex-direction: column;
    }}
}}

@media print {{
    body {{
        background: white;
    }}

    .report-container {{
        margin: 0;
        box-shadow: none;
    }}

    .disease-section,
    .figure-card {{
        break-inside: avoid;
    }}
}}
</style>
</head>

<body>
<div class="report-container">

<header class="cover">
    <div class="cover-label">
        Automated Variant Interpretation and Reporting
    </div>

    <h1>
        AI Genomics Platform:
        Complete Genomic Interpretation Report
    </h1>

    <p class="cover-subtitle">
        Integrated SNV/indel analysis, disease evidence,
        AI-assisted interpretation, comparative annotation,
        visualization and CNV workflow reporting.
    </p>

    <div class="cover-meta">
        Generated: {html.escape(generated)}
        &nbsp; | &nbsp;
        Research and educational use only
    </div>
</header>

<main class="content">

<section>
    <h2 class="section-title">1. Executive Summary</h2>

    <p class="lead">
        The AI Genomics Platform is a modular bioinformatics
        framework developed to standardize genomic variants,
        integrate disease and annotation evidence, generate
        explainable interpretations, visualize findings and
        export human-readable reports. The platform was
        evaluated using disease-associated variants in CFTR,
        HBB, MECP2 and FBN1 and includes a separate CNV
        integration branch.
    </p>

    {make_metric_cards(variants)}
</section>

<section>
    <h2 class="section-title">
        2. Platform Architecture and Workflow
    </h2>

    <div class="workflow">
        <div class="workflow-step">VCF Input</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">Normalization</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">VEP / SnpEff</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">Evidence Integration</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">AI Interpretation</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">Dashboard and Report</div>
    </div>

    <div class="workflow" style="margin-top:18px">
        <div class="workflow-step">CNV BED / SV-VCF</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">AnnotSV</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">ClassifyCNV</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">ISV-CNV / ClinGen</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">CNV Report</div>
    </div>
</section>

<section>
    <h2 class="section-title">
        3. Final Verified Variant Summary
    </h2>

    {dataframe_html(
        variant_summary,
        empty_message=(
            "The final verified variant table was unavailable."
        ),
    )}
</section>

<section>
    <h2 class="section-title">
        4. Disease-Specific Integrated Interpretation
    </h2>

    {make_disease_sections(
        variants,
        interpretations,
    )}
</section>

<section>
    <h2 class="section-title">
        5. Annotation and Database Tool Contributions
    </h2>

    {make_tools_table()}
</section>

<section>
    <h2 class="section-title">
        6. VEP, SnpEff and Reference Comparison
    </h2>

    {make_comparison_section()}
</section>

<section>
    <h2 class="section-title">
        7. AI-Assisted Interpretation Outputs
    </h2>

    {dataframe_html(
        interpretations,
        empty_message=(
            "AI interpretation output was unavailable."
        ),
        max_rows=20,
    )}
</section>

<section>
    <h2 class="section-title">
        8. Visual Analytics
    </h2>

    {make_figures_section()}
</section>

<section>
    <h2 class="section-title">
        9. CNV Workflow Results
    </h2>

    <div class="status-note">
        <strong>Implementation status:</strong>
        the CNV branch validates input, calculates CNV size,
        detects tool availability and generates CSV, JSON and
        HTML reports. AnnotSV, ClassifyCNV, ISV-CNV and ClinGen
        require external installation, resources or manual
        access before actual evidence generation.
    </div>

    <h3>CNV interpretation table</h3>

    {dataframe_html(
        cnvs,
        empty_message="CNV interpretation results unavailable.",
        max_rows=20,
    )}

    <h3>CNV tool status</h3>

    {dataframe_html(
        cnv_status,
        empty_message="CNV tool-status table unavailable.",
        max_rows=20,
    )}
</section>

<section>
    <h2 class="section-title">
        10. Strengths, Limitations and Future Development
    </h2>

    <h3>Key strengths</h3>
    <ul>
        <li>
            Modular SNV/indel and CNV workflow architecture.
        </li>
        <li>
            Integration of annotation, disease knowledge,
            visualization and AI-assisted explanation.
        </li>
        <li>
            Transparent handling of unavailable external data.
        </li>
        <li>
            Reproducible CSV, TSV, JSON and HTML outputs.
        </li>
        <li>
            Four disease-focused validation examples.
        </li>
    </ul>

    <h3>Current limitations</h3>
    <ul>
        <li>
            Some external database values were unavailable
            during the final offline run.
        </li>
        <li>
            A local VEP cache is required for repeat offline
            annotation.
        </li>
        <li>
            OMIM values are local reference metadata rather
            than live licensed API responses.
        </li>
        <li>
            AnnotSV, ClassifyCNV, ISV-CNV and ClinGen were not
            executed and remain external integrations.
        </li>
        <li>
            Automated outputs require professional genetics
            review before clinical use.
        </li>
    </ul>

    <h3>Future work</h3>
    <ul>
        <li>
            Restore and validate the local Ensembl VEP cache.
        </li>
        <li>
            Connect authenticated and licensed database APIs.
        </li>
        <li>
            Implement phenotype-driven Phen2Gene ranking.
        </li>
        <li>
            Install and execute the complete CNV tool chain.
        </li>
        <li>
            Expand ACMG evidence automation and provenance.
        </li>
    </ul>
</section>

<section>
    <h2 class="section-title">
        11. Conclusion
    </h2>

    <p class="lead">
        The AI Genomics Platform demonstrates an integrated,
        explainable and extensible approach to genomic variant
        interpretation. It combines variant standardization,
        annotation comparison, disease knowledge, AI-assisted
        reporting, visualization and a separate CNV framework.
        The project prioritizes transparency by clearly
        distinguishing available evidence from external
        resources that require further configuration.
    </p>
</section>

<section class="warning">
    <strong>Research-use warning:</strong>
    This report is intended for research and educational use.
    It is not a standalone diagnostic report. All findings must
    be reviewed by a qualified clinical genetics professional.
</section>

</main>

<footer class="footer">
    AI Genomics Platform —
    Automated Variant Interpretation and Reporting Framework
</footer>

</div>
</body>
</html>
"""


def main() -> None:
    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    SUBMISSION_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    report = build_report()

    OUTPUT_FILE.write_text(
        report,
        encoding="utf-8"    submission_copy = (
        SUBMISSION_DIR
        / OUTPUT_FILE.name
    )

    submission_copy.write_text(
        report,
        encoding="utf-8",
    )

    print("\nComplete master report generated successfully.")
    print(f"Main report: {OUTPUT_FILE}")
    print(f"Submission copy: {submission_copy}")


if __name__ == "__main__":
    main()
k
