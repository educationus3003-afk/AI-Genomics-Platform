#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate vep_env

mkdir -p ~/AI_Genomics_Platform/06_Results/VEP

vep \
  --offline \
  --cache \
  --species homo_sapiens \
  --cache_version 116 \
  --dir_cache /mnt/d/VEP_cache \
  --assembly GRCh38 \
  --input_file "$1" \
  --output_file "$2" \
  --vcf \
  --symbol \
  --canonical \
  --variant_class \
  --force_overwrite \
  --stats_file "${2%.vcf}_stats.html"
