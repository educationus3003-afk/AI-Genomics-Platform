from __future__ import annotations

import base64
import html
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


BASE = Path.home() / "AI_Genomics_Platform"

MASTER_TABLE = (
    BASE
    / "06_Results"
    / "Master"
    / "database_integrated_table.csv"
)

LIVE_TABLE = (
    BASE
    / "06_Results"
    / "Master"
    / "live_database_results.csv"
)

AI_TABLE = (
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

CNV_STATUS_TABLE = (
    BASE
    / "06_Results"
    / "CNV"
    / "cnv_tool_status.csv"
)

OUTPUT_REPORT = (
    BASE
    / "08_Reports"
    / "AI_Genomics_Platform_Complete_Report.html"
)

SUBMISSION_DIR = (
    BASE
    / "09_Presentation"
    / "Final_Submission"
)

SUBMISSION_REPORT = (
    SUBMISSION_DIR
    / "AI_Genomics_Platform_Complete_Report.html"
)

FIGURE_CANDIDATES = {
    "Gene Distribution": [
        BASE / "07_Visualizations" / "gene_distribution.png",
        BASE / "07_Visualizations" / "gene_distribution.jpg",
    ],
    "Database Coverage": [
        BASE / "07_Visualizations" / "database_coverage.png",
        BASE / "07_Visualizations" / "database_coverage.jpg",
    ],
    "Variant Classification": [
        BASE / "07_Visualizations" / "classification_pie.png",
        BASE / "07_Visualizations" / "classification_distribution.png",
    ],
}

MISSING_TEXT = {
    "",
    "nan",
    "none",
    "null",
    "not retrieved",
    "not available",
    "unavailable",
    "n/a",
}

DISEASE_INFORMATION = {
    "CFTR": {
        "disease": "Cystic Fibrosis",
        "inheritance": "Autosomal recessive",
        "overview": (
            "Cystic fibrosis is caused by disease-associated variants "
            "in CFTR. The gene encodes a chloride and bicarbonate channel "
            "that regulates epithelial ion and water transport."
        ),
        "mechanism": (
            "Reduced or absent CFTR activity produces dehydrated, thick "
            "secretions affecting the lungs, pancreas, gastrointestinal "
            "tract and reproductive system."
        ),
    },
    "HBB": {
        "disease": "Sickle Cell Disease / Beta-Globin Disorder",
        "inheritance": "Autosomal recessive",
        "overview": (
            "HBB encodes the beta-globin component of haemoglobin. "
            "Pathogenic variants can alter haemoglobin production or "
            "structure and cause inherited haemoglobin disorders."
        ),
        "mechanism": (
            "Abnormal beta-globin can cause haemoglobin polymerisation, "
            "red-cell deformation, haemolysis, vaso-occlusion or reduced "
            "beta-globin synthesis."
        ),
    },
    "MECP2": {
        "disease": "Rett Syndrome",
        "inheritance": "X-linked dominant",
        "overview": (
            "Rett syndrome is a neurodevelopmental disorder most commonly "
            "associated with pathogenic variants in MECP2."
        ),
        "mechanism": (
            "MECP2 dysfunction disrupts methylated-DNA binding, gene "
            "regulation, neuronal maturation and synaptic function."
        ),
    },
    "FBN1": {
        "disease": "Marfan Syndrome",
        "inheritance": "Autosomal dominant",
        "overview": (
            "Marfan syndrome is a connective-tissue disorder associated "
            "with pathogenic variants in FBN1, which encodes fibrillin-1."
        ),
        "mechanism": (
            "Abnormal fibrillin-1 affects extracellular microfibrils, "
            "tissue elasticity and regulation of TGF-beta signalling."
        ),
    },
}


def safe_read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()

    try:
        return pd.read_csv(path)
    except Exception as error:
        print(f"Warning: unable to read {path}: {error}")
        return pd.DataFrame()


def clean(value: Any, default: str = "") -> str:
    if value is None:
        return default

    try:
        if pd.isna(value):
            return default
    except (TypeError, ValueError):
        pass

    text = str(value).strip()

    if text.lower() in MISSING_TEXT:
        return default

    return text


def first_present(
    row: pd.Series,
    columns: Iterable[str],
    default: str = "Unavailable",
) -> str:
    for column in columns:
        if column not in row.index:
            continue

        value = clean(row[column])

        if value:
            return value

    return default


def standardise_variant_id(value: Any) -> str:
    text = clean(value)

    if not text:
        return ""

    text = text.replace(" ", "")
    return text


def load_combined_table() -> pd.DataFrame:
    master = safe_read_csv(MASTER_TABLE)
    live = safe_read_csv(LIVE_TABLE)

    if master.empty:
        raise FileNotFoundError(
            f"Master table was not found or was empty: {MASTER_TABLE}"
        )

    if "VARIANT_ID" not in master.columns:
        required = {"CHROM", "POS", "REF", "ALT"}

        if required.issubset(master.columns):
            master["VARIANT_ID"] = (
                "chr"
                + master["CHROM"]
                .astype(str)
                .str.replace("chr", "", regex=False)
                + ":"
                + master["POS"].astype(str)
                + ":"
                + master["REF"].astype(str)
                + ">"
                + master["ALT"].astype(str)
            )
        else:
            raise ValueError(
                "Master table requires VARIANT_ID or "
                "CHROM, POS, REF and ALT columns."
            )

    master["VARIANT_ID"] = master["VARIANT_ID"].map(
        standardise_variant_id
    )

    if "GENE" in master.columns:
        master["GENE"] = (
            master["GENE"]
            .astype(str)
            .str.upper()
            .str.strip()
        )

    if not live.empty and "VARIANT_ID" in live.columns:
        live["VARIANT_ID"] = live["VARIANT_ID"].map(
            standardise_variant_id
        )

        duplicate_columns = [
            column
            for column in live.columns
            if column in master.columns
            and column not in {"VARIANT_ID", "GENE"}
        ]

        master = master.drop(
            columns=duplicate_columns,
            errors="ignore",
        )

        merge_columns = [
            column
            for column in live.columns
            if column != "GENE"
        ]

        master = master.merge(
            live[merge_columns],
            on="VARIANT_ID",
            how="left",
        )

    return master


def image_data_uri(path: Path) -> str | None:
    if not path.exists():
        return None

    suffix = path.suffix.lower()

    mime_map = {
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".webp": "image/webp",
    }

    mime = mime_map.get(suffix)

    if mime is None:
        return None

    encoded = base64.b64encode(
        path.read_bytes()
    ).decode("ascii")

    return f"data:{mime};base64,{encoded}"


def table_html(
    dataframe: pd.DataFrame,
    empty_message: str,
    max_rows: int | None = None,
) -> str:
    if dataframe.empty:
        return (
            '<div class="empty-box">'
            f"{html.escape(empty_message)}"
            "</div>"
        )

    view = dataframe.copy()

    if max_rows is not None:
        view = view.head(max_rows)

    for column in view.columns:
        view[column] = view[column].map(
            lambda value: clean(value, "Unavailable")
        )

    return view.to_html(
        index=False,
        border=0,
        escape=True,
        classes="data-table",
    )


def database_status(value: str, database: str) -> str:
    cleaned = clean(value)

    if cleaned:
        return cleaned

    if database == "AlphaMissense":
        return "Not applicable or no eligible missense prediction returned"

    return f"No exact {database} result returned"


def get_ai_interpretation(
    gene: str,
    interpretations: pd.DataFrame,
) -> str:
    default = (
        "The evidence shown above should be evaluated together with "
        "phenotype, inheritance, population frequency, functional "
        "prediction, segregation and clinical database evidence."
    )

    if interpretations.empty:
        return default

    gene_column = next(
        (
            column
            for column in interpretations.columns
            if column.upper() == "GENE"
        ),
        None,
    )

    if gene_column is None:
        return default

    matching = interpretations[
        interpretations[gene_column]
        .astype(str)
        .str.upper()
        .str.strip()
        == gene
    ]

    if matching.empty:
        return default

    row = matching.iloc[0]

    return first_present(
        row,
        [
            "AI_Interpretation",
            "Interpretation",
            "Clinical_Interpretation",
            "AI_Summary",
            "Summary",
        ],
        default,
    )


def evidence_row(label: str, value: str) -> str:
    return (
        "<tr>"
        f"<th>{html.escape(label)}</th>"
        f"<td>{html.escape(value)}</td>"
        "</tr>"
    )


def build_disease_sections(
    variants: pd.DataFrame,
    interpretations: pd.DataFrame,
) -> str:
    output: list[str] = []

    for gene, details in DISEASE_INFORMATION.items():
        if "GENE" in variants.columns:
            rows = variants[
                variants["GENE"]
                .astype(str)
                .str.upper()
                .str.strip()
                == gene
            ]
        else:
            rows = pd.DataFrame()

        if rows.empty:
            output.append(
                f"""
                <section class="disease-card">
                    <div class="disease-header">
                        <div>
                            <span class="gene-badge">{gene}</span>
                            <h3>{html.escape(details["disease"])}</h3>
                        </div>
                    </div>

                    <div class="warning-box">
                        No matching variant record was found for {gene}
                        in the integrated master table.
                    </div>
                </section>
                """
            )
            continue

        for _, row in rows.iterrows():
            disease = first_present(
                row,
                [
                    "OMIM_DiseaseName",
                    "Disease",
                    "DISEASE",
                ],
                details["disease"],
            )

            inheritance = first_present(
                row,
                [
                    "OMIM_Inheritance",
                    "Inheritance",
                    "INHERITANCE",
                ],
                details["inheritance"],
            )

            variant_id = first_present(
                row,
                ["VARIANT_ID"],
            )

            consequence = first_present(
                row,
                [
                    "Live_VEP_Consequence",
                    "Final_Consequence",
                    "VEP_Consequence",
                    "Input_Consequence",
                    "CONSEQUENCE",
                ],
            )

            hgvsc = first_present(
                row,
                [
                    "Live_VEP_HGVSc",
                    "Final_HGVSc",
                    "VEP_HGVSc",
                    "Input_HGVSc",
                    "HGVSC",
                ],
            )

            hgvsp = first_present(
                row,
                [
                    "Live_VEP_HGVSp",
                    "Final_HGVSp",
                    "VEP_HGVSp",
                    "Input_HGVSp",
                    "HGVSP",
                ],
            )

            dbsnp = first_present(
                row,
                [
                    "Live_dbSNP_ID",
                    "Final_dbSNP_ID",
                    "dbSNP_ID",
                    "Input_dbSNP_ID",
                    "RSID",
                ],
                "",
            )

            clinvar = first_present(
                row,
                [
                    "Live_ClinVar",
                    "ClinVar_ClinicalSignificance",
                    "ClinVar",
                    "Input_ClinVar",
                ],
                "",
            )

            clinvar_review = first_present(
                row,
                [
                    "Live_ClinVar_ReviewStatus",
                    "ClinVar_ReviewStatus",
                ],
                "",
            )

            gnomad = first_present(
                row,
                [
                    "Live_gnomAD_AF",
                    "gnomAD_AlleleFrequency",
                    "gnomAD_AF",
                    "Input_gnomAD_AF",
                ],
                "",
            )

            alpha_prediction = first_present(
                row,
                [
                    "Live_AlphaMissense_Prediction",
                    "AlphaMissense_Prediction",
                    "Input_AlphaMissense",
                ],
                "",
            )

            alpha_score = first_present(
                row,
                [
                    "Live_AlphaMissense_Score",
                    "AlphaMissense_Score",
                ],
                "",
            )

            omim_id = first_present(
                row,
                [
                    "OMIM_OMIM_ID",
                    "OMIM_ID",
                ],
            )

            evidence_source = first_present(
                row,
                [
                    "Evidence_Sources",
                    "Annotation_Source",
                    "VEP_Source",
                    "SOURCE",
                ],
                "Integrated project evidence",
            )

            ensembl_status = first_present(
                row,
                ["Ensembl_Status"],
                "Status not recorded",
            )

            myvariant_status = first_present(
                row,
                ["MyVariant_Status"],
                "Status not recorded",
            )

            alpha_display = database_status(
                alpha_prediction,
                "AlphaMissense",
            )

            if alpha_score:
                alpha_display = (
                    f"{alpha_display}; score: {alpha_score}"
                )

            evidence = "".join(
                [
                    evidence_row("Variant", variant_id),
                    evidence_row("Gene", gene),
                    evidence_row("Disease", disease),
                    evidence_row("Inheritance", inheritance),
                    evidence_row("Consequence", consequence),
                    evidence_row("HGVSc", hgvsc),
                    evidence_row("HGVSp", hgvsp),
                    evidence_row(
                        "dbSNP",
                        database_status(dbsnp, "dbSNP"),
                    ),
                    evidence_row(
                        "ClinVar significance",
                        database_status(
                            clinvar,
                            "ClinVar",
                        ),
                    ),
                    evidence_row(
                        "ClinVar review status",
                        database_status(
                            clinvar_review,
                            "ClinVar review-status",
                        ),
                    ),
                    evidence_row(
                        "gnomAD allele frequency",
                        database_status(
                            gnomad,
                            "gnomAD frequency",
                        ),
                    ),
                    evidence_row(
                        "AlphaMissense",
                        alpha_display,
                    ),
                    evidence_row("OMIM ID", omim_id),
                    evidence_row(
                        "Evidence sources",
                        evidence_source,
                    ),
                    evidence_row(
                        "Ensembl query status",
                        ensembl_status,
                    ),
                    evidence_row(
                        "MyVariant query status",
                        myvariant_status,
                    ),
                ]
            )

            interpretation = get_ai_interpretation(
                gene,
                interpretations,
            )

            output.append(
                f"""
                <section class="disease-card">
                    <div class="disease-header">
                        <div>
                            <span class="gene-badge">{gene}</span>
                            <h3>{html.escape(disease)}</h3>
                        </div>

                        <span class="inheritance-badge">
                            {html.escape(inheritance)}
                        </span>
                    </div>

                    <div class="description-grid">
                        <div class="description-box">
                            <h4>Disease overview</h4>
                            <p>{html.escape(details["overview"])}</p>
                        </div>

                        <div class="description-box">
                            <h4>Molecular mechanism</h4>
                            <p>{html.escape(details["mechanism"])}</p>
                        </div>
                    </div>

                    <h4>Integrated variant evidence</h4>

                    <table class="evidence-table">
                        {evidence}
                    </table>

                    <div class="interpretation-box">
                        <h4>AI-assisted interpretation</h4>
                        <p>{html.escape(interpretation)}</p>
                    </div>
                </section>
                """
            )

    return "\n".join(output)


def build_figure_section() -> str:
    cards: list[str] = []

    for title, paths in FIGURE_CANDIDATES.items():
        selected = next(
            (path for path in paths if path.exists()),
            None,
        )

        if selected is None:
            continue

        uri = image_data_uri(selected)

        if uri is None:
            continue

        cards.append(
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

    if not cards:
        return (
            '<div class="empty-box">'
            "No visualization image files were available."
            "</div>"
        )

    return (
        '<div class="figure-grid">'
        + "".join(cards)
        + "</div>"
    )


def build_coverage_table(variants: pd.DataFrame) -> pd.DataFrame:
    mappings = {
        "VEP consequence": [
            "Live_VEP_Consequence",
            "Final_Consequence",
            "VEP_Consequence",
        ],
        "HGVS coding": [
            "Live_VEP_HGVSc",
            "Final_HGVSc",
            "VEP_HGVSc",
        ],
        "HGVS protein": [
            "Live_VEP_HGVSp",
            "Final_HGVSp",
            "VEP_HGVSp",
        ],
        "dbSNP": [
            "Live_dbSNP_ID",
            "Final_dbSNP_ID",
            "dbSNP_ID",
        ],
        "ClinVar": [
            "Live_ClinVar",
            "ClinVar_ClinicalSignificance",
        ],
        "gnomAD": [
            "Live_gnomAD_AF",
            "gnomAD_AF",
        ],
        "AlphaMissense": [
            "Live_AlphaMissense_Prediction",
            "AlphaMissense_Prediction",
        ],
        "OMIM": [
            "OMIM_OMIM_ID",
            "OMIM_ID",
        ],
    }

    records = []

    for database, columns in mappings.items():
        available = 0

        for _, row in variants.iterrows():
            value = first_present(
                row,
                columns,
                "",
            )

            if clean(value):
                available += 1

        records.append(
            {
                "Evidence type": database,
                "Available records": available,
                "Total variants": len(variants),
                "Coverage": (
                    f"{(available / len(variants) * 100):.1f}%"
                    if len(variants)
                    else "0.0%"
                ),
            }
        )

    return pd.DataFrame(records)


def build_report() -> str:
    variants = load_combined_table()
    interpretations = safe_read_csv(AI_TABLE)
    cnv_results = safe_read_csv(CNV_TABLE)
    cnv_status = safe_read_csv(CNV_STATUS_TABLE)

    generated = datetime.now().strftime(
        "%d %B %Y at %H:%M"
    )

    preferred_columns = [
        "VARIANT_ID",
        "GENE",
        "Live_dbSNP_ID",
        "Live_ClinVar",
        "Live_ClinVar_ReviewStatus",
        "Live_gnomAD_AF",
        "Live_VEP_Consequence",
        "Live_VEP_HGVSc",
        "Live_VEP_HGVSp",
        "Live_AlphaMissense_Prediction",
        "Live_AlphaMissense_Score",
        "OMIM_OMIM_ID",
        "OMIM_DiseaseName",
        "OMIM_Inheritance",
        "Evidence_Sources",
    ]

    summary_columns = [
        column
        for column in preferred_columns
        if column in variants.columns
    ]

    summary = (
        variants[summary_columns]
        if summary_columns
        else variants
    )

    coverage = build_coverage_table(variants)

    gene_count = (
        variants["GENE"].nunique()
        if "GENE" in variants.columns
        else 0
    )

    disease_sections = build_disease_sections(
        variants,
        interpretations,
    )

    figures = build_figure_section()

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta
    name="viewport"
    content="width=device-width, initial-scale=1.0"
>
<title>AI Genomics Platform Complete Report</title>

<style>
:root {{
    --navy: #102a43;
    --blue: #1769aa;
    --light-blue: #eaf4fb;
    --green: #18794e;
    --light-green: #edfdf5;
    --orange: #b45309;
    --light-orange: #fff7ed;
    --red: #b42318;
    --light-red: #fff1f0;
    --grey: #52606d;
    --light-grey: #f4f7fa;
    --border: #d7e0e8;
    --white: #ffffff;
}}

* {{
    box-sizing: border-box;
}}

body {{
    margin: 0;
    background: #edf2f7;
    color: #243b53;
    font-family: "Segoe UI", Arial, sans-serif;
    line-height: 1.6;
}}

.report {{
    max-width: 1450px;
    margin: 25px auto;
    background: var(--white);
    box-shadow: 0 12px 38px rgba(16, 42, 67, 0.14);
}}

.cover {{
    padding: 72px 70px 60px;
    color: white;
    background:
        linear-gradient(
            135deg,
            #102a43,
            #1769aa 70%,
            #3490c7
        );
}}

.cover-label {{
    display: inline-block;
    padding: 7px 15px;
    border: 1px solid rgba(255,255,255,0.45);
    border-radius: 20px;
    font-size: 13px;
    font-weight: 700;
    letter-spacing: 1px;
}}

.cover h1 {{
    max-width: 1000px;
    margin: 22px 0 15px;
    font-size: 46px;
    line-height: 1.15;
}}

.cover p {{
    max-width: 1000px;
    margin: 0;
    color: #e1eff8;
    font-size: 19px;
}}

.cover-meta {{
    margin-top: 32px;
    color: #d4e8f5;
    font-size: 14px;
}}

.main-content {{
    padding: 45px 65px 70px;
}}

.main-section {{
    margin-bottom: 55px;
}}

.section-title {{
    margin: 0 0 22px;
    padding-bottom: 10px;
    color: var(--navy);
    border-bottom: 3px solid var(--blue);
    font-size: 29px;
}}

.lead {{
    max-width: 1150px;
    color: #425466;
    font-size: 17px;
}}

.metrics {{
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 18px;
    margin-top: 28px;
}}

.metric {{
    padding: 24px;
    border: 1px solid var(--border);
    border-radius: 12px;
    background: var(--light-grey);
}}

.metric-number {{
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
    gap: 9px;
    margin-top: 24px;
}}

.workflow-step {{
    padding: 12px 17px;
    border-radius: 8px;
    background: var(--light-blue);
    color: var(--navy);
    font-weight: 700;
}}

.workflow-arrow {{
    color: var(--blue);
    font-size: 23px;
    font-weight: 800;
}}

.data-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 13px;
    overflow-wrap: anywhere;
}}

.data-table th {{
    padding: 11px 10px;
    background: var(--navy);
    color: white;
    text-align: left;
    vertical-align: top;
}}

.data-table td {{
    padding: 10px;
    border-bottom: 1px solid var(--border);
    vertical-align: top;
}}

.data-table tr:nth-child(even) {{
    background: #f7f9fb;
}}

.disease-card {{
    margin-bottom: 30px;
    padding: 27px;
    border: 1px solid var(--border);
    border-radius: 14px;
    box-shadow: 0 5px 18px rgba(16,42,67,0.07);
}}

.disease-header {{
    display: flex;
    justify-content: space-between;
    align-items: flex-start;
    gap: 18px;
    margin-bottom: 22px;
}}

.disease-header h3 {{
    margin: 8px 0 0;
    color: var(--navy);
    font-size: 25px;
}}

.gene-badge {{
    display: inline-block;
    padding: 5px 12px;
    border-radius: 6px;
    background: var(--blue);
    color: white;
    font-weight: 800;
}}

.inheritance-badge {{
    padding: 7px 13px;
    border-radius: 18px;
    background: var(--light-green);
    color: var(--green);
    font-size: 13px;
    font-weight: 700;
}}

.description-grid {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 17px;
    margin-bottom: 24px;
}}

.description-box {{
    padding: 17px;
    border-left: 4px solid var(--blue);
    background: var(--light-grey);
}}

.description-box h4 {{
    margin: 0 0 7px;
    color: var(--navy);
}}

.description-box p {{
    margin: 0;
}}

.evidence-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 14px;
}}

.evidence-table th {{
    width: 250px;
    padding: 10px;
    background: #f0f4f8;
    color: var(--navy);
    border: 1px solid var(--border);
    text-align: left;
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
    background: var(--light-green);
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

.notice-box {{
    padding: 18px;
    border-left: 5px solid var(--orange);
    background: var(--light-orange);
}}

.warning-box {{
    padding: 18px;
    border-left: 5px solid var(--red);
    background: var(--light-red);
}}

.empty-box {{
    padding: 16px;
    border: 1px dashed var(--border);
    background: var(--light-grey);
    color: var(--grey);
}}

.footer {{
    padding: 25px 65px;
    background: var(--navy);
    color: #d7e7f3;
    font-size: 13px;
}}

@media (max-width: 900px) {{
    .cover {{
        padding: 45px 27px;
    }}

    .cover h1 {{
        font-size: 34px;
    }}

    .main-content {{
        padding: 35px 23px;
    }}

    .metrics,
    .description-grid,
    .figure-grid {{
        grid-template-columns: 1fr;
    }}

    .disease-header {{
        flex-direction: column;
    }}
}}

@media print {{
    body {{
        background: white;
    }}

    .report {{
        margin: 0;
        box-shadow: none;
    }}

    .disease-card,
    .figure-card {{
        break-inside: avoid;
    }}
}}
</style>
</head>

<body>
<div class="report">

<header class="cover">
    <span class="cover-label">
        Integrated Variant Interpretation Platform
    </span>

    <h1>
        AI Genomics Platform:
        Complete Genomic Interpretation Report
    </h1>

    <p>
        A consolidated report combining variant annotation,
        database evidence, disease interpretation,
        AI-assisted explanations, visual analytics and
        copy-number variant workflow results.
    </p>

    <div class="cover-meta">
        Generated on {html.escape(generated)}
        &nbsp; | &nbsp;
        Research and educational use only
    </div>
</header>

<main class="main-content">

<section class="main-section">
    <h2 class="section-title">1. Executive Summary</h2>

    <p class="lead">
        The AI Genomics Platform integrates genomic variant
        standardization, functional annotation, disease knowledge,
        population evidence, clinical database evidence,
        AI-assisted interpretation, visualization and reporting.
        The current validation set includes CFTR, HBB, MECP2 and
        FBN1 disease-associated variants.
    </p>

    <div class="metrics">
        <div class="metric">
            <div class="metric-number">{len(variants)}</div>
            <div class="metric-label">Integrated variants</div>
        </div>

        <div class="metric">
            <div class="metric-number">{gene_count}</div>
            <div class="metric-label">Genes evaluated</div>
        </div>

        <div class="metric">
            <div class="metric-number">4</div>
            <div class="metric-label">Disease models</div>
        </div>

        <div class="metric">
            <div class="metric-number">2</div>
            <div class="metric-label">SNV/indel and CNV branches</div>
        </div>
    </div>
</section>

<section class="main-section">
    <h2 class="section-title">
        2. Platform Workflow
    </h2>

    <div class="workflow">
        <div class="workflow-step">VCF Input</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">Normalization</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">VEP / SnpEff</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">Database Integration</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">AI Interpretation</div>
        <div class="workflow-arrow">→</div>
        <div class="workflow-step">Complete HTML Report</div>
    </div>

    <div class="workflow">
        <div class="workflow-step">CNV Input</div>
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

<section class="main-section">
    <h2 class="section-title">
        3. Final Integrated Variant Table
    </h2>

    {table_html(
        summary,
        "No integrated variant records were available.",
    )}
</section>

<section class="main-section">
    <h2 class="section-title">
        4. Database Evidence Coverage
    </h2>

    {table_html(
        coverage,
        "Database coverage could not be calculated.",
    )}
</section>

<section class="main-section">
    <h2 class="section-title">
        5. Disease-Specific Interpretation
    </h2>

    {disease_sections}
</section>

<section class="main-section">
    <h2 class="section-title">
        6. AI Interpretation Output Table
    </h2>

    {table_html(
        interpretations,
        "AI interpretation output was unavailable.",
        max_rows=30,
    )}
</section>

<section class="main-section">
    <h2 class="section-title">
        7. Visual Analytics
    </h2>

    {figures}
</section>

<section class="main-section">
    <h2 class="section-title">
        8. Copy-Number Variant Branch
    </h2>

    <div class="notice-box">
        The CNV branch supports input validation, CNV-size
        calculation, tool-status detection and report generation.
        AnnotSV, ClassifyCNV, ISV-CNV and ClinGen evidence require
        their corresponding external installations or resources.
    </div>

    <h3>CNV interpretation results</h3>

    {table_html(
        cnv_results,
        "No CNV interpretation records were available.",
        max_rows=30,
    )}

    <h3>CNV tool status</h3>

    {table_html(
        cnv_status,
        "No CNV tool-status records were available.",
        max_rows=30,
    )}
</section>

<section class="main-section">
    <h2 class="section-title">
        9. Limitations and Future Improvements
    </h2>

    <ul>
        <li>
            A missing exact ClinVar or gnomAD result does not
            indicate that a variant is benign.
        </li>
        <li>
            AlphaMissense is applicable only to eligible
            missense substitutions.
        </li>
        <li>
            A local Ensembl VEP cache is required for reliable
            repeat offline annotation.
        </li>
        <li>
            Live online database queries require an active
            internet connection.
        </li>
        <li>
            OMIM licensed API access may be required for
            automated live OMIM retrieval.
        </li>
        <li>
            Clinical interpretation requires professional
            review and phenotype correlation.
        </li>
    </ul>
</section>

<section class="main-section">
    <h2 class="section-title">10. Conclusion</h2>

    <p class="lead">
        The AI Genomics Platform provides a unified and
        transparent framework for genomic variant interpretation.
        It consolidates available functional, clinical,
        population and disease evidence into one professional
        report while clearly distinguishing retrieved evidence,
        unavailable matches and predictions that are not
        scientifically applicable.
    </p>
</section>

<section class="warning-box">
    <strong>Research-use warning:</strong>
    This report is intended for research and educational use.
    It is not a standalone clinical diagnosis. All findings
    require review by a qualified genetics professional.
</section>

</main>

<footer class="footer">
    AI Genomics Platform —
    Automated Genomic Variant Interpretation Framework
</footer>

</div>
</body>
</html>
"""


def main() -> None:
    OUTPUT_REPORT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    SUBMISSION_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    report = build_report()

    OUTPUT_REPORT.write_text(
        report,
        encoding="utf-8",
    )

    SUBMISSION_REPORT.write_text(
        report,
        encoding="utf-8",
    )

    print("\nComplete report generated successfully.")
    print(f"Main report: {OUTPUT_REPORT}")
    print(f"Submission copy: {SUBMISSION_REPORT}")


if __name__ == "__main__":
    main()
