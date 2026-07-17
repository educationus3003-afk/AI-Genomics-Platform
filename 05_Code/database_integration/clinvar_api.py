from pathlib import Path
import pandas as pd

BASE = Path.home() / "AI_Genomics_Platform"

INPUT = BASE / "06_Results/Master/database_integrated_table.csv"
OUTPUT = BASE / "06_Results/Master/database_integrated_table.csv"


# Demo mapping until live API is connected
CLINVAR = {
    "CFTR": {
        "ClinVar_ClinicalSignificance": "Pathogenic",
        "ClinVar_ReviewStatus": "Practice guideline",
        "ClinVar_Disease": "Cystic Fibrosis",
    },

    "HBB": {
        "ClinVar_ClinicalSignificance": "Pathogenic",
        "ClinVar_ReviewStatus": "Practice guideline",
        "ClinVar_Disease": "Beta Thalassemia",
    },

    "MECP2": {
        "ClinVar_ClinicalSignificance": "Likely pathogenic",
        "ClinVar_ReviewStatus": "Expert panel",
        "ClinVar_Disease": "Rett Syndrome",
    },

    "FBN1": {
        "ClinVar_ClinicalSignificance": "Pathogenic",
        "ClinVar_ReviewStatus": "Expert panel",
        "ClinVar_Disease": "Marfan Syndrome",
    }
}


df = pd.read_csv(INPUT)

for i, row in df.iterrows():

    gene = row["GENE"]

    if gene in CLINVAR:

        for column, value in CLINVAR[gene].items():
            df.loc[i, column] = value

df.to_csv(OUTPUT, index=False)

print(df[
[
"GENE",
"ClinVar_ClinicalSignificance",
"ClinVar_ReviewStatus",
"ClinVar_Disease"
]
])
