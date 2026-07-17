from __future__ import annotations

import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


BASE = Path.home() / "AI_Genomics_Platform"
LOG_DIR = BASE / "06_Results" / "Logs"
LOG_FILE = LOG_DIR / "ai_genomics_pipeline.log"

STEPS = [
    (
        "Standardize annotations",
        BASE / "05_Code/annotation_engine/annotation_engine.py",
    ),
    (
        "Create database evidence table",
        BASE / "05_Code/database_integration/database_integrator.py",
    ),
    (
        "Extract real VEP evidence",
        BASE / "05_Code/database_integration/extract_vep_evidence.py",
    ),
    (
        "Enrich variants with dbSNP, gnomAD and AlphaMissense",
        BASE / "05_Code/database_integration/variant_enrichment.py",
    ),
    (
        "Generate AI-assisted interpretations",
        BASE / "05_Code/ai_engine/ai_interpreter.py",
    ),
    (
        "Generate visualization dashboard",
        BASE / "05_Code/visualization/dashboard_generator.py",
    ),
    (
        "Generate HTML report",
        BASE / "05_Code/report_generator/report_generator.py",
    ),
]

REQUIRED_INPUTS = [
    BASE / "06_Results/Master/master_annotation_table.csv",
    BASE / "06_Results/VEP/sample_input_vep.vcf",
    BASE / "01_Disease_Knowledge/diseases.json",
]


def write_log(message: str) -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {message}"

    print(line)

    with LOG_FILE.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def validate_inputs() -> None:
    missing = [path for path in REQUIRED_INPUTS if not path.exists()]

    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(
            "Required pipeline input files are missing:\n"
            f"{formatted}"
        )


def run_step(number: int, title: str, script: Path) -> None:
    if not script.exists():
        raise FileNotFoundError(f"Pipeline module not found: {script}")

    write_log(f"STEP {number}/{len(STEPS)} STARTED: {title}")

    start = time.perf_counter()

    result = subprocess.run(
        [sys.executable, str(script)],
        cwd=BASE,
        capture_output=True,
        text=True,
        check=False,
    )

    elapsed = time.perf_counter() - start

    if result.stdout.strip():
        with LOG_FILE.open("a", encoding="utf-8") as handle:
            handle.write(result.stdout.rstrip() + "\n")

    if result.stderr.strip():
        with LOG_FILE.open("a", encoding="utf-8") as handle:
            handle.write("STDERR:\n")
            handle.write(result.stderr.rstrip() + "\n")

    if result.returncode != 0:
        write_log(
            f"STEP {number} FAILED: {title} "
            f"(exit code {result.returncode})"
        )

        print("\nLast error output:\n")
        print(result.stderr or result.stdout)

        raise RuntimeError(f"Pipeline stopped during: {title}")

    write_log(
        f"STEP {number} COMPLETED: {title} "
        f"({elapsed:.2f} seconds)"
    )


def display_outputs() -> None:
    outputs = [
        BASE / "06_Results/Master/database_integrated_table.csv",
        BASE / "06_Results/Interpretation/ai_interpretations.json",
        BASE / "06_Results/Interpretation/ai_interpretations.csv",
        BASE / "07_Visualizations/interactive_dashboard.html",
        BASE / "07_Visualizations/gene_distribution.png",
        BASE / "07_Visualizations/classification_pie.png",
        BASE / "07_Visualizations/database_coverage.png",
        BASE / "08_Reports/genomic_ai_report.html",
    ]

    print("\nGenerated outputs:\n")

    for path in outputs:
        status = "READY" if path.exists() else "MISSING"
        print(f"[{status}] {path}")


def main() -> None:
    pipeline_start = time.perf_counter()

    print("\n" + "=" * 72)
    print("AI GENOMICS PLATFORM")
    print("Automated Variant Interpretation and Reporting Pipeline")
    print("=" * 72 + "\n")

    LOG_DIR.mkdir(parents=True, exist_ok=True)

    with LOG_FILE.open("a", encoding="utf-8") as handle:
        handle.write("\n" + "=" * 72 + "\n")

    write_log("Pipeline execution started")

    try:
        validate_inputs()

        for number, (title, script) in enumerate(STEPS, start=1):
            run_step(number, title, script)

    except Exception as error:
        write_log(f"PIPELINE FAILED: {error}")
        sys.exit(1)

    total_time = time.perf_counter() - pipeline_start

    write_log(
        f"Pipeline completed successfully in {total_time:.2f} seconds"
    )

    display_outputs()

    print("\n" + "=" * 72)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("=" * 72)
    print(f"\nLog file: {LOG_FILE}")
    print(
        "\nResearch-use warning: automated results require review by "
        "a qualified genetics professional."
    )


if __name__ == "__main__":
    main()
