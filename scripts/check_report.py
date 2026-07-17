from pathlib import Path
import hashlib
import sys


MASTER = Path(
    "08_Reports/AI_Genomics_Platform_Complete_Report.html"
)

PAGES = Path("docs/index.html")

required_terms = [
    "AI Genomics Platform",
    "Executive Summary",
    "Platform Workflow",
    "CFTR",
    "HBB",
    "MECP2",
    "FBN1",
    "ClinVar",
    "gnomAD",
    "AI-assisted interpretation",
    "Research",
]

forbidden_terms = [
    "GeneInsightAI",
    "AI Genomic Platform",
]

failed = False


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


for path in (MASTER, PAGES):
    if not path.exists():
        print(f"ERROR: Missing file: {path}")
        failed = True

if failed:
    sys.exit(1)

master_text = MASTER.read_text(
    encoding="utf-8",
    errors="replace",
)

for term in required_terms:
    if term.lower() not in master_text.lower():
        print(f"MISSING REQUIRED CONTENT: {term}")
        failed = True

for term in forbidden_terms:
    if term.lower() in master_text.lower():
        print(f"FOUND OLD PROJECT NAME: {term}")
        failed = True

if digest(MASTER) != digest(PAGES):
    print("ERROR: Master report and docs/index.html are different.")
    failed = True

if "not validated for clinical" not in master_text.lower():
    print(
        "WARNING: Add an explicit statement that the platform "
        "is not validated for clinical use."
    )
    failed = True

if failed:
    sys.exit(1)

print("Report quality checks passed.")
print(f"Master: {MASTER}")
print(f"Pages:  {PAGES}")
print(f"SHA256: {digest(MASTER)}")
