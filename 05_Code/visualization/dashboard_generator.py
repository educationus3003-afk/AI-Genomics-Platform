from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px

BASE = Path.home() / "AI_Genomics_Platform"

INPUT = BASE / "06_Results" / "Master" / "database_integrated_table.csv"
OUTPUT = BASE / "07_Visualizations"

OUTPUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

# -------------------------
# Gene Distribution
# -------------------------
gene_counts = df["GENE"].value_counts()

plt.figure(figsize=(6,4))
gene_counts.plot(kind="bar")
plt.title("Gene Distribution")
plt.ylabel("Variant Count")
plt.tight_layout()
plt.savefig(OUTPUT/"gene_distribution.png")
plt.close()

# -------------------------
# Classification Distribution
# -------------------------
classification = df["INTERVAR_CLASSIFICATION"].fillna("Unknown").value_counts()

plt.figure(figsize=(6,6))
classification.plot(kind="pie", autopct="%1.1f%%")
plt.ylabel("")
plt.title("Variant Classification")
plt.tight_layout()
plt.savefig(OUTPUT/"classification_pie.png")
plt.close()

# -------------------------
# Database Coverage
# -------------------------
columns = [
    "ClinVar_ClinicalSignificance",
    "gnomAD_AlleleFrequency",
    "AlphaMissense_Prediction",
    "REVEL_Score",
    "CADD_PHRED",
    "SpliceAI_DeltaScore"
]

coverage = {}

for c in columns:
    if c in df.columns:
        coverage[c] = (df[c] != "Not Retrieved").sum()

coverage_df = pd.DataFrame(
    {
        "Database": list(coverage.keys()),
        "Retrieved": list(coverage.values())
    }
)

plt.figure(figsize=(8,4))
plt.bar(coverage_df["Database"], coverage_df["Retrieved"])
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(OUTPUT/"database_coverage.png")
plt.close()

# -------------------------
# Interactive Plotly Dashboard
# -------------------------
fig = px.sunburst(
    df,
    path=["GENE","INTERVAR_CLASSIFICATION"],
    title="Genomic Variant Dashboard"
)

fig.write_html(OUTPUT/"interactive_dashboard.html")

print("\nDashboard created successfully.\n")
print(OUTPUT)
