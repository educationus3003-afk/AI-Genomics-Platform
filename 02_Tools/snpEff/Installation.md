Tool: SnpEff

Purpose:
Predict the functional effects of genomic variants and annotate each variant according to its biological consequence.

Installation:

conda install -c bioconda snpeff -y

Genome Database:

GRCh38.99

Commands Used:

snpEff --version

snpEff databases

snpEff GRCh38.99 sample_input.vcf > sample_input_snpeff.vcf

Output:

Annotated VCF containing ANN field.

Status:

SUCCESS
