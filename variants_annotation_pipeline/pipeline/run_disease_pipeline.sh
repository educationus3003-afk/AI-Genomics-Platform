#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'

# Find the project root automatically
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

usage() {
    echo "Usage:"
    echo "  bash pipeline/run_disease_pipeline.sh <disease_name>"
    echo
    echo "Examples:"
    echo "  bash pipeline/run_disease_pipeline.sh fabry"
    echo "  bash pipeline/run_disease_pipeline.sh monilethrix"
    echo "  bash pipeline/run_disease_pipeline.sh ais"
    echo "  bash pipeline/run_disease_pipeline.sh nf1"
}

# Disease name is the first command-line argument
DISEASE="${1:-}"

if [[ -z "$DISEASE" ]]; then
    usage
    exit 1
fi

# Permit safe filename characters only
if [[ ! "$DISEASE" =~ ^[a-z0-9_]+$ ]]; then
    echo "ERROR: Use lowercase letters, numbers or underscores only."
    exit 1
fi

SNV_INPUT="$PROJECT_ROOT/input/small_variants/${DISEASE}.small_variants.vcf"
CNV_INPUT="$PROJECT_ROOT/input/cnv/${DISEASE}.cnvs.bed"
OUTDIR="$PROJECT_ROOT/results/${DISEASE}"

if [[ ! -s "$SNV_INPUT" ]]; then
    echo "ERROR: Missing SNV input:"
    echo "$SNV_INPUT"
    exit 1
fi

if [[ ! -s "$CNV_INPUT" ]]; then
    echo "ERROR: Missing CNV input:"
    echo "$CNV_INPUT"
    exit 1
fi

mkdir -p \
    "$OUTDIR/snv" \
    "$OUTDIR/cnv" \
    "$OUTDIR/work" \
    "$OUTDIR/logs" \
    "$OUTDIR/final"

echo "Input validation successful"
echo "Disease: $DISEASE"
echo "SNV input: $SNV_INPUT"
echo "CNV input: $CNV_INPUT"
echo "Output directory: $OUTDIR"

# ============================================================
# Container and reference paths
# ============================================================
CORE_SIF="$PROJECT_ROOT/containers/core_tools.sif"
REF_FASTA="$PROJECT_ROOT/resources/reference/hg38.fa"

[[ -s "$CORE_SIF" ]] || {
    echo "ERROR: Missing container: $CORE_SIF"
    exit 1
}

[[ -s "$REF_FASTA" ]] || {
    echo "ERROR: Missing reference: $REF_FASTA"
    exit 1
}

[[ -s "${REF_FASTA}.fai" ]] || {
    echo "ERROR: Missing reference index: ${REF_FASTA}.fai"
    exit 1
}

# Paths as seen inside the container
CONTAINER_SNV_INPUT="/project/input/small_variants/${DISEASE}.small_variants.vcf"
CONTAINER_OUTDIR="/project/results/${DISEASE}"

NORMALIZED_VCF="$OUTDIR/work/${DISEASE}.normalized.vcf.gz"
CONTAINER_NORMALIZED_VCF="${CONTAINER_OUTDIR}/work/${DISEASE}.normalized.vcf.gz"

LOGFILE="$OUTDIR/logs/${DISEASE}.pipeline.log"

exec > >(tee -a "$LOGFILE") 2>&1

echo
echo "============================================================"
echo "STEP 1: SNV/INDEL NORMALIZATION"
echo "============================================================"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$CORE_SIF" \
    bcftools norm \
    -f /project/resources/reference/hg38.fa \
    -m -any \
    -c x \
    -Oz \
    -o "$CONTAINER_NORMALIZED_VCF" \
    "$CONTAINER_SNV_INPUT"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$CORE_SIF" \
    tabix -f -p vcf "$CONTAINER_NORMALIZED_VCF"

NORMALIZED_RECORDS=$(
    apptainer exec \
        --bind "$PROJECT_ROOT:/project" \
        "$CORE_SIF" \
        bcftools view -H "$CONTAINER_NORMALIZED_VCF" |
    wc -l
)

if [[ "$NORMALIZED_RECORDS" -eq 0 ]]; then
    echo "ERROR: Normalization produced an empty VCF."
    exit 1
fi

echo "PASS: Normalization completed"
echo "Variant records retained: $NORMALIZED_RECORDS"
echo "Output: $NORMALIZED_VCF"

# ============================================================
# STEP 2: VEP ANNOTATION
# ============================================================
VEP_SIF="$PROJECT_ROOT/containers/vep.sif"
VEP_CACHE="$PROJECT_ROOT/resources/vep_cache"

[[ -s "$VEP_SIF" ]] || {
    echo "ERROR: Missing VEP container: $VEP_SIF"
    exit 1
}

[[ -d "$VEP_CACHE/homo_sapiens/115_GRCh38" ]] || {
    echo "ERROR: VEP 115 GRCh38 cache not found."
    echo "Expected: $VEP_CACHE/homo_sapiens/115_GRCh38"
    exit 1
}

VEP_OUTPUT="$OUTDIR/snv/${DISEASE}.vep.annotated.vcf"
CONTAINER_VEP_OUTPUT="${CONTAINER_OUTDIR}/snv/${DISEASE}.vep.annotated.vcf"

echo
echo "============================================================"
echo "STEP 2: VEP 115 ANNOTATION"
echo "============================================================"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$VEP_SIF" \
    vep \
    --input_file "$CONTAINER_NORMALIZED_VCF" \
    --output_file "$CONTAINER_VEP_OUTPUT" \
    --format vcf \
    --vcf \
    --offline \
    --cache \
    --dir_cache /project/resources/vep_cache \
    --species homo_sapiens \
    --assembly GRCh38 \
    --fasta /project/resources/reference/hg38.fa \
    --symbol \
    --canonical \
    --mane \
    --hgvs \
    --variant_class \
    --biotype \
    --numbers \
    --force_overwrite \
    --no_stats \
    --fork 4

if [[ ! -s "$VEP_OUTPUT" ]]; then
    echo "ERROR: VEP output was not created."
    exit 1
fi

if ! grep -q '^##INFO=<ID=CSQ' "$VEP_OUTPUT"; then
    echo "ERROR: VEP CSQ annotation is missing."
    exit 1
fi

VEP_RECORDS=$(
    apptainer exec \
        --bind "$PROJECT_ROOT:/project" \
        "$CORE_SIF" \
        bcftools view -H "$CONTAINER_VEP_OUTPUT" |
    wc -l
)

echo "PASS: VEP annotation completed"
echo "Variant records: $VEP_RECORDS"
echo "Annotation field present: CSQ"
echo "Output: $VEP_OUTPUT"

# ============================================================
# STEP 2: VEP ANNOTATION
# ============================================================
VEP_SIF="$PROJECT_ROOT/containers/vep.sif"
VEP_CACHE="$PROJECT_ROOT/resources/vep_cache"

[[ -s "$VEP_SIF" ]] || {
    echo "ERROR: Missing VEP container: $VEP_SIF"
    exit 1
}

[[ -d "$VEP_CACHE/homo_sapiens/115_GRCh38" ]] || {
    echo "ERROR: VEP 115 GRCh38 cache not found."
    echo "Expected: $VEP_CACHE/homo_sapiens/115_GRCh38"
    exit 1
}

VEP_OUTPUT="$OUTDIR/snv/${DISEASE}.vep.annotated.vcf"
CONTAINER_VEP_OUTPUT="${CONTAINER_OUTDIR}/snv/${DISEASE}.vep.annotated.vcf"

echo
echo "============================================================"
echo "STEP 2: VEP 115 ANNOTATION"
echo "============================================================"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$VEP_SIF" \
    vep \
    --input_file "$CONTAINER_NORMALIZED_VCF" \
    --output_file "$CONTAINER_VEP_OUTPUT" \
    --format vcf \
    --vcf \
    --offline \
    --cache \
    --dir_cache /project/resources/vep_cache \
    --species homo_sapiens \
    --assembly GRCh38 \
    --fasta /project/resources/reference/hg38.fa \
    --symbol \
    --canonical \
    --mane \
    --hgvs \
    --variant_class \
    --biotype \
    --numbers \
    --force_overwrite \
    --no_stats \
    --fork 4

if [[ ! -s "$VEP_OUTPUT" ]]; then
    echo "ERROR: VEP output was not created."
    exit 1
fi

if ! grep -q '^##INFO=<ID=CSQ' "$VEP_OUTPUT"; then
    echo "ERROR: VEP CSQ annotation is missing."
    exit 1
fi

VEP_RECORDS=$(
    apptainer exec \
        --bind "$PROJECT_ROOT:/project" \
        "$CORE_SIF" \
        bcftools view -H "$CONTAINER_VEP_OUTPUT" |
    wc -l
)

echo "PASS: VEP annotation completed"
echo "Variant records: $VEP_RECORDS"
echo "Annotation field present: CSQ"
echo "Output: $VEP_OUTPUT"

# ============================================================
# STEP 3: SNPEFF ANNOTATION
# ============================================================
SNPEFF_SIF="$PROJECT_ROOT/containers/snpeff.sif"
SNPEFF_DATA="$PROJECT_ROOT/resources/snpeff_data/data"
SNPEFF_GENOME="GRCh38.mane.1.2.ensembl"

[[ -s "$SNPEFF_SIF" ]] || {
    echo "ERROR: Missing SnpEff container: $SNPEFF_SIF"
    exit 1
}

[[ -d "$SNPEFF_DATA/$SNPEFF_GENOME" ]] || {
    echo "ERROR: SnpEff database not found:"
    echo "$SNPEFF_DATA/$SNPEFF_GENOME"
    exit 1
}

SNPEFF_OUTPUT="$OUTDIR/snv/${DISEASE}.vep.snpeff.annotated.vcf.gz"
CONTAINER_SNPEFF_OUTPUT="${CONTAINER_OUTDIR}/snv/${DISEASE}.vep.snpeff.annotated.vcf.gz"

echo
echo "============================================================"
echo "STEP 3: SNPEFF ANNOTATION"
echo "============================================================"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$SNPEFF_SIF" \
    /opt/snpEff/exec/snpeff ann \
    -noStats \
    -canon \
    -hgvs \
    -dataDir /project/resources/snpeff_data/data \
    "$SNPEFF_GENOME" \
    "$CONTAINER_VEP_OUTPUT" \
| apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$CORE_SIF" \
    bgzip -c \
> "$SNPEFF_OUTPUT"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$CORE_SIF" \
    tabix -f -p vcf "$CONTAINER_SNPEFF_OUTPUT"

ANN_HEADER=$(
    apptainer exec \
        --bind "$PROJECT_ROOT:/project" \
        "$CORE_SIF" \
        bcftools view -h "$CONTAINER_SNPEFF_OUTPUT" |
    grep -m1 '^##INFO=<ID=ANN' || true
)

if [[ -z "$ANN_HEADER" ]]; then
    echo "ERROR: SnpEff ANN annotation is missing."
    exit 1
fi

SNPEFF_RECORDS=$(
    apptainer exec \
        --bind "$PROJECT_ROOT:/project" \
        "$CORE_SIF" \
        bcftools view -H "$CONTAINER_SNPEFF_OUTPUT" |
    wc -l
)

echo "PASS: SnpEff annotation completed"
echo "Variant records: $SNPEFF_RECORDS"
echo "Annotation field present: ANN"
echo "Output: $SNPEFF_OUTPUT"

# ============================================================
# STEP 4: CLINVAR ANNOTATION
# ============================================================
CLINVAR_VCF="$PROJECT_ROOT/resources/clinvar/clinvar.chr.vcf.gz"
CLINVAR_INDEX="${CLINVAR_VCF}.tbi"

[[ -s "$CLINVAR_VCF" ]] || {
    echo "ERROR: Missing ClinVar database:"
    echo "$CLINVAR_VCF"
    exit 1
}

[[ -s "$CLINVAR_INDEX" ]] || {
    echo "ERROR: Missing ClinVar index:"
    echo "$CLINVAR_INDEX"
    exit 1
}

CLINVAR_OUTPUT="$OUTDIR/snv/${DISEASE}.vep.snpeff.clinvar.annotated.vcf.gz"
CONTAINER_CLINVAR_OUTPUT="${CONTAINER_OUTDIR}/snv/${DISEASE}.vep.snpeff.clinvar.annotated.vcf.gz"

echo
echo "============================================================"
echo "STEP 4: CLINVAR ANNOTATION"
echo "============================================================"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$CORE_SIF" \
    bcftools annotate \
    -a /project/resources/clinvar/clinvar.chr.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNDISDB,INFO/CLNHGVS,INFO/CLNVC,INFO/CLNVCSO,INFO/GENEINFO \
    -Oz \
    -o "$CONTAINER_CLINVAR_OUTPUT" \
    "$CONTAINER_SNPEFF_OUTPUT"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$CORE_SIF" \
    tabix -f -p vcf "$CONTAINER_CLINVAR_OUTPUT"

if [[ ! -s "$CLINVAR_OUTPUT" ]]; then
    echo "ERROR: ClinVar output was not created."
    exit 1
fi

if ! apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$CORE_SIF" \
    bcftools view -h "$CONTAINER_CLINVAR_OUTPUT" |
    grep -q '^##INFO=<ID=CLNSIG'
then
    echo "ERROR: ClinVar CLNSIG header is missing."
    exit 1
fi

CLINVAR_MATCHES=$(
    apptainer exec \
        --bind "$PROJECT_ROOT:/project" \
        "$CORE_SIF" \
        bcftools query \
        -f '%INFO/CLNSIG\n' \
        "$CONTAINER_CLINVAR_OUTPUT" |
    awk '$1 != "." && $1 != "" {count++} END {print count+0}'
)

echo "PASS: ClinVar annotation completed"
echo "Variants with ClinVar CLNSIG matches: $CLINVAR_MATCHES"
echo "Output: $CLINVAR_OUTPUT"

# ============================================================
# STEP 5: SPLICEAI ANNOTATION
# ============================================================
SPLICEAI_SIF="$PROJECT_ROOT/containers/spliceai.sif"

[[ -s "$SPLICEAI_SIF" ]] || {
    echo "ERROR: Missing SpliceAI container:"
    echo "$SPLICEAI_SIF"
    exit 1
}

SPLICEAI_RAW="$OUTDIR/work/${DISEASE}.spliceai.raw.vcf"
SPLICEAI_OUTPUT="$OUTDIR/snv/${DISEASE}.vep.snpeff.clinvar.spliceai.vcf.gz"

CONTAINER_SPLICEAI_RAW="${CONTAINER_OUTDIR}/work/${DISEASE}.spliceai.raw.vcf"
CONTAINER_SPLICEAI_OUTPUT="${CONTAINER_OUTDIR}/snv/${DISEASE}.vep.snpeff.clinvar.spliceai.vcf.gz"

echo
echo "============================================================"
echo "STEP 5: SPLICEAI ANNOTATION"
echo "============================================================"

rm -f \
    "$SPLICEAI_RAW" \
    "$SPLICEAI_OUTPUT" \
    "${SPLICEAI_OUTPUT}.tbi"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$SPLICEAI_SIF" \
    spliceai \
    -I "$CONTAINER_CLINVAR_OUTPUT" \
    -O "$CONTAINER_SPLICEAI_RAW" \
    -R /project/resources/reference/hg38.fa \
    -A grch38

if [[ ! -s "$SPLICEAI_RAW" ]]; then
    echo "ERROR: SpliceAI did not create an output VCF."
    exit 1
fi

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$CORE_SIF" \
    bgzip -c "$CONTAINER_SPLICEAI_RAW" \
    > "$SPLICEAI_OUTPUT"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$CORE_SIF" \
    tabix -f -p vcf "$CONTAINER_SPLICEAI_OUTPUT"

if ! apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$CORE_SIF" \
    bcftools view -h "$CONTAINER_SPLICEAI_OUTPUT" |
    grep -q '^##INFO=<ID=SpliceAI'
then
    echo "ERROR: SpliceAI annotation header is missing."
    exit 1
fi

SPLICEAI_MATCHES=$(
    apptainer exec \
        --bind "$PROJECT_ROOT:/project" \
        "$CORE_SIF" \
        bcftools query \
        -f '%INFO/SpliceAI\n' \
        "$CONTAINER_SPLICEAI_OUTPUT" |
    awk '$1 != "." && $1 != "" {count++} END {print count+0}'
)

echo "PASS: SpliceAI annotation completed"
echo "Variants with SpliceAI annotation: $SPLICEAI_MATCHES"
echo "Output: $SPLICEAI_OUTPUT"

# ============================================================
# STEP 6: ANNOTSV CNV ANNOTATION
# ============================================================
ANNOTSV_SIF="$PROJECT_ROOT/containers/annotsv.sif"
ANNOTSV_DATABASE="$PROJECT_ROOT/resources/annotsv_annotations/AnnotSV_annotations"

[[ -s "$ANNOTSV_SIF" ]] || {
    echo "ERROR: Missing AnnotSV container:"
    echo "$ANNOTSV_SIF"
    exit 1
}

[[ -d "$ANNOTSV_DATABASE/Annotations_Human" ]] || {
    echo "ERROR: AnnotSV human annotation database not found:"
    echo "$ANNOTSV_DATABASE/Annotations_Human"
    exit 1
}

ANNOTSV_OUTPUT="$OUTDIR/cnv/${DISEASE}.AnnotSV.tsv"
CONTAINER_ANNOTSV_OUTPUT="${CONTAINER_OUTDIR}/cnv/${DISEASE}.AnnotSV.tsv"
CONTAINER_CNV_INPUT="/project/input/cnv/${DISEASE}.cnvs.bed"

echo
echo "============================================================"
echo "STEP 6: ANNOTSV CNV ANNOTATION"
echo "============================================================"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$ANNOTSV_SIF" \
    AnnotSV \
    -SVinputFile "$CONTAINER_CNV_INPUT" \
    -outputFile "$CONTAINER_ANNOTSV_OUTPUT" \
    -annotationsDir /project/resources/annotsv_annotations/AnnotSV_annotations \
    -genomeBuild GRCh38 \
    -svtBEDcol 4 \
    -annotationMode both \
    -overwrite 1

if [[ ! -s "$ANNOTSV_OUTPUT" ]]; then
    echo "ERROR: AnnotSV output was not created."
    exit 1
fi

if ! head -n 1 "$ANNOTSV_OUTPUT" | grep -q 'AnnotSV_ID'; then
    echo "ERROR: AnnotSV output header is missing."
    exit 1
fi

ANNOTSV_FULL_ROWS=$(
    awk -F '\t' '
        NR==1 {
            for (i=1; i<=NF; i++)
                if ($i=="Annotation_mode") col=i
            next
        }
        col && $col=="full" {count++}
        END {print count+0}
    ' "$ANNOTSV_OUTPUT"
)

echo "PASS: AnnotSV annotation completed"
echo "Whole-CNV records annotated: $ANNOTSV_FULL_ROWS"
echo "Output: $ANNOTSV_OUTPUT"

# ============================================================
# STEP 7: CLASSIFYCNV ACMG/CLINGEN CLASSIFICATION
# ============================================================
CLASSIFYCNV_DIR="$PROJECT_ROOT/tools/ClassifyCNV"
CLASSIFYCNV_SCRIPT="$CLASSIFYCNV_DIR/ClassifyCNV.py"

[[ -s "$CLASSIFYCNV_SCRIPT" ]] || {
    echo "ERROR: ClassifyCNV.py not found:"
    echo "$CLASSIFYCNV_SCRIPT"
    exit 1
}

# ClassifyCNV previously worked when its output directory was created
# relative to the ClassifyCNV installation folder.
CLASSIFYCNV_RUN_NAME="${DISEASE}_ClassifyCNV"
CLASSIFYCNV_RUN_DIR="$CLASSIFYCNV_DIR/$CLASSIFYCNV_RUN_NAME"

CLASSIFYCNV_OUTPUT="$OUTDIR/cnv/${DISEASE}.ClassifyCNV.Scoresheet.txt"

echo
echo "============================================================"
echo "STEP 7: CLASSIFYCNV CLASSIFICATION"
echo "============================================================"

# Remove an older run so repeated pipeline executions remain clean
rm -rf "$CLASSIFYCNV_RUN_DIR"
rm -f "$CLASSIFYCNV_OUTPUT"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    --pwd /project/tools/ClassifyCNV \
    "$CORE_SIF" \
    python3 ClassifyCNV.py \
    --infile "$CONTAINER_CNV_INPUT" \
    --GenomeBuild hg38 \
    --cores 4 \
    --precise \
    --outdir "$CLASSIFYCNV_RUN_NAME"

CLASSIFYCNV_SCORESHEET="$CLASSIFYCNV_RUN_DIR/Scoresheet.txt"

if [[ ! -s "$CLASSIFYCNV_SCORESHEET" ]]; then
    echo "ERROR: ClassifyCNV Scoresheet.txt was not created."
    exit 1
fi

cp "$CLASSIFYCNV_SCORESHEET" "$CLASSIFYCNV_OUTPUT"

if ! head -n 1 "$CLASSIFYCNV_OUTPUT" | grep -q 'Classification'; then
    echo "ERROR: ClassifyCNV classification column is missing."
    exit 1
fi

CLASSIFYCNV_ROWS=$(
    awk 'NR > 1 && NF > 0 {count++} END {print count+0}' \
        "$CLASSIFYCNV_OUTPUT"
)

echo "PASS: ClassifyCNV completed"
echo "CNV records classified: $CLASSIFYCNV_ROWS"
echo "Output: $CLASSIFYCNV_OUTPUT"

# ============================================================
# STEP 8: ISV-CNV PREDICTION
# ============================================================
ISV_SIF="$PROJECT_ROOT/containers/isv.sif"

[[ -s "$ISV_SIF" ]] || {
    echo "ERROR: Missing ISV container:"
    echo "$ISV_SIF"
    exit 1
}

ISV_INPUT="$OUTDIR/work/${DISEASE}.isv.bed"
ISV_OUTPUT="$OUTDIR/cnv/${DISEASE}.ISV_with_SHAP.tsv"

CONTAINER_ISV_INPUT="${CONTAINER_OUTDIR}/work/${DISEASE}.isv.bed"
CONTAINER_ISV_OUTPUT="${CONTAINER_OUTDIR}/cnv/${DISEASE}.ISV_with_SHAP.tsv"

echo
echo "============================================================"
echo "STEP 8: ISV-CNV PREDICTION"
echo "============================================================"

# ISV requires a header row above the four CNV BED columns
{
    printf "chromosome\tstart\tend\tcnv_type\n"
    cat "$CNV_INPUT"
} > "$ISV_INPUT"

rm -f "$ISV_OUTPUT"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$ISV_SIF" \
    python3 - \
    "$CONTAINER_ISV_INPUT" \
    "$CONTAINER_ISV_OUTPUT" <<'PY'
import sys

import pandas as pd
import isv

input_file = sys.argv[1]
output_file = sys.argv[2]

cnvs = pd.read_csv(input_file, sep="\t")

required_columns = [
    "chromosome",
    "start",
    "end",
    "cnv_type"
]

if list(cnvs.columns) != required_columns:
    raise ValueError(
        f"Expected columns {required_columns}, "
        f"but found {list(cnvs.columns)}"
    )

if cnvs.empty:
    raise ValueError("ISV input contains no CNVs.")

results = isv.isv(
    cnvs,
    proba=True,
    shap=True,
    threshold=0.95
)

results.to_csv(
    output_file,
    sep="\t",
    index=False
)
PY

if [[ ! -s "$ISV_OUTPUT" ]]; then
    echo "ERROR: ISV-CNV output was not created."
    exit 1
fi

if ! head -n 1 "$ISV_OUTPUT" | grep -q $'\tISV\t'; then
    echo "ERROR: ISV probability column is missing."
    exit 1
fi

ISV_ROWS=$(
    awk 'NR > 1 && NF > 0 {count++} END {print count+0}' \
        "$ISV_OUTPUT"
)

echo "PASS: ISV-CNV completed"
echo "CNV records predicted: $ISV_ROWS"
echo "Output: $ISV_OUTPUT"

# ============================================================
# STEP 9: COMBINE CNV RESULTS
# ============================================================
CNV_SUMMARY="$OUTDIR/final/${DISEASE}.final.cnv.summary.tsv"
CONTAINER_CNV_SUMMARY="${CONTAINER_OUTDIR}/final/${DISEASE}.final.cnv.summary.tsv"

echo
echo "============================================================"
echo "STEP 9: COMBINING CNV RESULTS"
echo "============================================================"

rm -f "$CNV_SUMMARY"

apptainer exec \
    --bind "$PROJECT_ROOT:/project" \
    "$ISV_SIF" \
    python3 - \
    "$CONTAINER_ANNOTSV_OUTPUT" \
    "/project/results/${DISEASE}/cnv/${DISEASE}.ClassifyCNV.Scoresheet.txt" \
    "$CONTAINER_ISV_OUTPUT" \
    "$CONTAINER_CNV_SUMMARY" <<'PY'
import sys
import pandas as pd

annotsv_file = sys.argv[1]
classify_file = sys.argv[2]
isv_file = sys.argv[3]
output_file = sys.argv[4]

keys = ["Chromosome", "Start", "End", "Type"]

# ---------------- AnnotSV ----------------
ann = pd.read_csv(
    annotsv_file,
    sep="\t",
    dtype=str,
    low_memory=False
)

ann = ann[ann["Annotation_mode"] == "full"].copy()

ann["Chromosome"] = (
    "chr" + ann["SV_chrom"].str.replace("chr", "", regex=False)
)

# AnnotSV reports the BED start as one-based
ann["Start"] = pd.to_numeric(
    ann["SV_start"],
    errors="raise"
) - 1

ann["End"] = pd.to_numeric(
    ann["SV_end"],
    errors="raise"
)

ann["Type"] = ann["SV_type"]

acmg_labels = {
    1: "Benign",
    2: "Likely benign",
    3: "VUS",
    4: "Likely pathogenic",
    5: "Pathogenic"
}

ann["AnnotSV_ACMG_class"] = pd.to_numeric(
    ann["ACMG_class"],
    errors="coerce"
)

ann["AnnotSV_classification"] = (
    ann["AnnotSV_ACMG_class"].map(acmg_labels)
)

ann = ann[
    keys + [
        "Gene_count",
        "AnnotSV_ranking_score",
        "AnnotSV_ACMG_class",
        "AnnotSV_classification"
    ]
].drop_duplicates(subset=keys)

# ---------------- ClassifyCNV ----------------
classify = pd.read_csv(
    classify_file,
    sep="\t",
    dtype=str
)

classify["Start"] = pd.to_numeric(
    classify["Start"],
    errors="raise"
)

classify["End"] = pd.to_numeric(
    classify["End"],
    errors="raise"
)

classify = classify[
    keys + [
        "Classification",
        "Total score"
    ]
].rename(
    columns={
        "Classification": "ClassifyCNV_classification",
        "Total score": "ClassifyCNV_total_score"
    }
)

# ---------------- ISV-CNV ----------------
isv = pd.read_csv(
    isv_file,
    sep="\t"
)

isv = isv.rename(
    columns={
        "chrom": "Chromosome",
        "start": "Start",
        "end": "End",
        "cnv_type": "Type",
        "ISV": "ISV_probability"
    }
)

isv["ISV_threshold_0.95"] = isv["ISV_probability"].apply(
    lambda value:
        "Above threshold"
        if value >= 0.95
        else "Below threshold"
)

isv = isv[
    keys + [
        "ISV_probability",
        "ISV_threshold_0.95"
    ]
]

# ---------------- Merge ----------------
combined = (
    ann.merge(
        classify,
        on=keys,
        how="outer"
    )
    .merge(
        isv,
        on=keys,
        how="outer"
    )
    .sort_values(keys)
)

combined.to_csv(
    output_file,
    sep="\t",
    index=False
)

print(combined.to_string(index=False))
PY

if [[ ! -s "$CNV_SUMMARY" ]]; then
    echo "ERROR: Combined CNV summary was not created."
    exit 1
fi

CNV_SUMMARY_ROWS=$(
    awk 'NR > 1 && NF > 0 {count++} END {print count+0}' \
        "$CNV_SUMMARY"
)

echo "PASS: Combined CNV summary created"
echo "CNV records summarized: $CNV_SUMMARY_ROWS"
echo "Output: $CNV_SUMMARY"

# ============================================================
# STEP 10: CREATE FINAL SNV RESULT
# ============================================================
FINAL_SNV="$OUTDIR/final/${DISEASE}.final.small_variants.annotated.vcf.gz"
FINAL_SNV_INDEX="${FINAL_SNV}.tbi"

echo
echo "============================================================"
echo "STEP 10: CREATING FINAL SNV RESULT"
echo "============================================================"

rm -f "$FINAL_SNV" "$FINAL_SNV_INDEX"

cp "$SPLICEAI_OUTPUT" "$FINAL_SNV"

# The index is optional for submission, but useful for fast querying
if [[ -s "${SPLICEAI_OUTPUT}.tbi" ]]; then
    cp "${SPLICEAI_OUTPUT}.tbi" "$FINAL_SNV_INDEX"
fi

if [[ ! -s "$FINAL_SNV" ]]; then
    echo "ERROR: Final SNV file was not created."
    exit 1
fi

FINAL_SNV_RECORDS=$(
    apptainer exec \
        --bind "$PROJECT_ROOT:/project" \
        "$CORE_SIF" \
        bcftools view -H \
        "/project/results/${DISEASE}/final/${DISEASE}.final.small_variants.annotated.vcf.gz" |
    wc -l
)

FINAL_SNV_HEADER=$(
    apptainer exec \
        --bind "$PROJECT_ROOT:/project" \
        "$CORE_SIF" \
        bcftools view -h \
        "/project/results/${DISEASE}/final/${DISEASE}.final.small_variants.annotated.vcf.gz"
)

for annotation in CSQ ANN CLNSIG SpliceAI
do
    if grep -q "^##INFO=<ID=${annotation}" <<< "$FINAL_SNV_HEADER"; then
        echo "PASS annotation: $annotation"
    else
        echo "WARNING: Annotation header not found: $annotation"
    fi
done

echo "PASS: Final SNV result created"
echo "Variant records: $FINAL_SNV_RECORDS"
echo "Output: $FINAL_SNV"

if [[ -s "$FINAL_SNV_INDEX" ]]; then
    echo "Optional index: $FINAL_SNV_INDEX"
fi

# ============================================================
# STEP 11: FINAL VALIDATION AND OUTPUT MANIFEST
# ============================================================
FINAL_MANIFEST="$OUTDIR/final/${DISEASE}.pipeline_outputs.txt"

echo
echo "============================================================"
echo "STEP 11: FINAL PIPELINE VALIDATION"
echo "============================================================"

REQUIRED_OUTPUTS=(
    "$FINAL_SNV"
    "$CNV_SUMMARY"
    "$ANNOTSV_OUTPUT"
    "$CLASSIFYCNV_OUTPUT"
    "$ISV_OUTPUT"
)

VALIDATION_FAILED=0

for file in "${REQUIRED_OUTPUTS[@]}"
do
    if [[ -s "$file" ]]; then
        echo "PASS: $file"
    else
        echo "ERROR: Missing or empty output: $file"
        VALIDATION_FAILED=1
    fi
done

if [[ "$VALIDATION_FAILED" -ne 0 ]]; then
    echo "ERROR: Final pipeline validation failed."
    exit 1
fi

{
    echo "RARE-DISEASE ANNOTATION PIPELINE OUTPUTS"
    echo "========================================"
    echo
    echo "Disease: $DISEASE"
    echo "Genome build: GRCh38"
    echo "Completion time: $(date --iso-8601=seconds)"
    echo
    echo "INPUT FILES"
    echo "SNV: $SNV_INPUT"
    echo "CNV: $CNV_INPUT"
    echo
    echo "FINAL OUTPUTS"
    echo "Annotated SNV VCF: $FINAL_SNV"

    if [[ -s "$FINAL_SNV_INDEX" ]]; then
        echo "Optional SNV index: $FINAL_SNV_INDEX"
    fi

    echo "Combined CNV summary: $CNV_SUMMARY"
    echo
    echo "DETAILED CNV OUTPUTS"
    echo "AnnotSV: $ANNOTSV_OUTPUT"
    echo "ClassifyCNV: $CLASSIFYCNV_OUTPUT"
    echo "ISV-CNV: $ISV_OUTPUT"
    echo
    echo "ANNOTATION TOOLS"
    echo "- bcftools normalization"
    echo "- VEP 115"
    echo "- SnpEff"
    echo "- ClinVar"
    echo "- SpliceAI"
    echo "- AnnotSV"
    echo "- ClassifyCNV"
    echo "- ISV-CNV"
    echo
    echo "SNV records: $FINAL_SNV_RECORDS"
    echo "CNV records: $CNV_SUMMARY_ROWS"
    echo
    echo "Pipeline log: $LOGFILE"
} > "$FINAL_MANIFEST"

echo
echo "============================================================"
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "============================================================"
echo "Disease: $DISEASE"
echo "Final SNV: $FINAL_SNV"
echo "Final CNV: $CNV_SUMMARY"
echo "Manifest: $FINAL_MANIFEST"
echo "Log: $LOGFILE"
