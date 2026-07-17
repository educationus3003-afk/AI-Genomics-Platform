# GPT vs Pipeline Annotation Comparison

**Diseases:** Fabry, Monilethrix, Androgen Insensitivity Syndrome (AIS), and NF1  
**Updated:** AIS pipeline `VCF.GZ` included; the previous tabix-index-only entry was replaced.

## Final Verdict

- **16/16 usable data files**
- **8/8 GPT-versus-pipeline comparisons completed**
- **8/8 core disease-causing findings agree**
- **4/8 whole-file comparisons are exact matches**
  - All four CNV coordinate sets match.
  - The complete SNV record sets do not match.

The main biological conclusions agree for all four diseases. The same disease-causing SNVs and CNVs are identified, with consistent main gene, consequence, HGVS, ClinVar, AnnotSV, and ClassifyCNV interpretations where applicable.

However, each pipeline SNV VCF contains **6 variants**, while each GPT SNV VCF contains **19 variants**. The 6 retained variants include the disease-causing variant and match the GPT files in coordinates and basic VCF fields, but **13 GPT background variants are absent** from every pipeline SNV output.

---

## Pair-Level Summary

| Disease | Core SNV finding | Whole SNV file | Core CNV finding | Whole CNV coordinate set | Interpretation |
|---|---|---|---|---|---|
| Fabry | MATCH | MISMATCH | MATCH | MATCH | The disease-causing SNV and CNV agree, but the complete SNV set contains only 6 of 19 GPT records. |
| Monilethrix | MATCH | MISMATCH | MATCH | MATCH | The disease-causing SNV and CNV agree, but the complete SNV set contains only 6 of 19 GPT records. |
| AIS | MATCH | MISMATCH | MATCH | MATCH | The AR SNV and deletion agree after replacing the tabix index with the actual pipeline VCF.GZ. |
| NF1 | MATCH | MISMATCH | MATCH | MATCH | The disease-causing SNV and CNV agree, but the complete SNV set contains only 6 of 19 GPT records. |

---

## Uploaded-File Inventory

### GPT-Annotated Files

| Type | Disease | File |
|---|---|---|
| SNV | Fabry | `SYN_FABRY.final.small_variants.annotated.vcf` |
| CNV | Fabry | `SYN_FABRY.final.cnv.annotated.vcf` |
| SNV | Monilethrix | `SYN_MONILETHRIX.final.small_variants.annotated.vcf` |
| CNV | Monilethrix | `SYN_MONILETHRIX.final.cnv.annotated.vcf` |
| SNV | AIS | `SYN_AIS.final.small_variants.annotated.vcf` |
| CNV | AIS | `SYN_AIS.final.cnv.annotated(1).vcf` |
| SNV | NF1 | `SYN_NF1.final.small_variants.annotated.vcf` |
| CNV | NF1 | `SYN_NF1.final.cnv.annotated.vcf` |

### Pipeline-Annotated Files

| Type | Disease | File |
|---|---|---|
| SNV | Fabry | `fabry.final.small_variants.annotated.vcf.gz` |
| CNV | Fabry | `fabry.final.cnv.summary.tsv` |
| SNV | Monilethrix | `monilethrix.final.small_variants.annotated.vcf.gz` |
| CNV | Monilethrix | `monilethrix.final.cnv.summary.tsv` |
| SNV | AIS | `ais.final.small_variants.annotated.vcf.gz` |
| CNV | AIS | `ais.final.cnv.summary.tsv` |
| SNV | NF1 | `nf1.final.small_variants.annotated.vcf.gz` |
| CNV | NF1 | `nf1.final.cnv.summary.tsv` |

---

## Small-Variant Record-Set Comparison

| Disease | GPT records | Pipeline records | Exact overlap | GPT records missing from pipeline | Pipeline-only records | Whole-file status |
|---|---:|---:|---:|---:|---:|---|
| Fabry | 19 | 6 | 6 | 13 | 0 | MISMATCH |
| Monilethrix | 19 | 6 | 6 | 13 | 0 | MISMATCH |
| AIS | 19 | 6 | 6 | 13 | 0 | MISMATCH |
| NF1 | 19 | 6 | 6 | 13 | 0 | MISMATCH |

For every disease:

- The 6 shared records match in chromosome, position, reference allele, alternate allele, ID, QUAL, FILTER, FORMAT, and genotype.
- The pipeline output contains no additional variants that are absent from the GPT file.
- The 13 missing records appear to be synthetic background variants.
- The exact reason for their removal should be confirmed from the pipeline normalization or filtering logs before stating a definite cause.

---

## Disease-Causing SNV Comparison

### Fabry Disease

- **Variant:** `chrX:101398467 C>T`
- **Gene:** `GLA`
- **Protein change:** `p.Arg301Gln`
- **Result:** GPT and pipeline annotations agree.
- **Status:** **MATCH**

### Monilethrix

- **Variant:** `chr12:52306272 G>T`
- **Gene:** `KRT86`
- **Protein change:** `p.Glu413Asp`
- **Result:** GPT and pipeline annotations agree.
- **Status:** **MATCH**

### Androgen Insensitivity Syndrome

- **Variant:** `chrX:67721873 C>T`
- **Gene:** `AR`
- **Protein change:** `p.Arg787Ter`
- **Result:** GPT and pipeline annotations agree after the actual pipeline VCF.GZ replaced the tabix index.
- **Status:** **MATCH**

### Neurofibromatosis Type 1

- **Variant:** `chr17:31230383 G>A`
- **Gene:** `NF1`
- **Coding effect:** `c.3113+1G>A`
- **Consequence:** splice-donor variant
- **Pipeline evidence:** SpliceAI supports donor loss.
- **Result:** GPT and pipeline annotations agree.
- **Status:** **MATCH**

The GPT files contain compact curated annotation text, whereas the pipeline files contain fuller transcript-level VEP, SnpEff, ClinVar, and SpliceAI annotations. Therefore, INFO fields are not expected to be identical character-for-character.

---

## CNV Coordinate-Set Comparison

| Disease | GPT CNVs | Pipeline CNVs | Coordinate/type matches | Status |
|---|---:|---:|---:|---|
| Fabry | 9 | 9 | 9 | MATCH |
| Monilethrix | 9 | 9 | 9 | MATCH |
| AIS | 9 | 9 | 9 | MATCH |
| NF1 | 9 | 9 | 9 | MATCH |

All four CNV pairs contain the same 9 CNVs.

The apparent one-base difference at the start of each CNV interval is expected:

- VCF `POS` is **1-based**.
- The pipeline TSV `Start` is **BED-style 0-based**.

Therefore, a VCF start of `1001` corresponds to a TSV start of `1000`.

---

## Disease-Causing CNV Comparison

### Fabry Disease

- **Target gene:** `GLA`
- **Result:** The GPT and pipeline deletion coordinates and interpretations agree.
- **Status:** **MATCH**

### Monilethrix

- **Target gene:** `DSG4`
- **Result:** The deletion coordinates agree.
- **Tool-level difference:** AnnotSV reports **Pathogenic**, while ClassifyCNV reports **Uncertain significance**.
- This is not a GPT-versus-pipeline mismatch because the GPT file records the same disagreement between tools.
- **Status:** **MATCH for file comparison; tool interpretations remain discordant**

### Androgen Insensitivity Syndrome

- **Target gene:** `AR`
- **Result:** The GPT and pipeline deletion coordinates and interpretations agree.
- **Status:** **MATCH**

### Neurofibromatosis Type 1

- **Target gene:** `NF1`
- **Result:** The GPT and pipeline deletion coordinates and interpretations agree.
- **Pipeline evidence:** ISV-CNV exceeds the reported pathogenicity threshold.
- **Status:** **MATCH**

---

## Disease-Wise Conclusions

### Fabry

The `GLA` p.Arg301Gln SNV and the `GLA` deletion agree between GPT and pipeline.

### Monilethrix

The `KRT86` p.Glu413Asp SNV agrees. The `DSG4` deletion coordinates also agree, including the same AnnotSV-versus-ClassifyCNV disagreement.

### AIS

The `AR` p.Arg787Ter SNV now matches successfully, and the `AR` deletion also agrees between GPT and pipeline.

### NF1

The `NF1` c.3113+1G>A splice-donor variant agrees, including the pipeline SpliceAI evidence. The `NF1` deletion also agrees.

---

## Overall Conclusion

The comparison demonstrates that the **main disease-causing findings are consistent between GPT annotation and the automated pipeline for all four diseases**.

The CNV results match completely after accounting for the expected 1-based versus 0-based coordinate convention.

The remaining limitation is that every pipeline SNV file contains only **6 of the 19 GPT variants**. Therefore:

- **Core biological interpretation:** MATCH
- **CNV files:** MATCH
- **Disease-causing SNVs:** MATCH
- **Complete SNV record sets:** MISMATCH

The missing background variants should be investigated using the pipeline normalization and filtering logs before the final report explains why they were removed.
