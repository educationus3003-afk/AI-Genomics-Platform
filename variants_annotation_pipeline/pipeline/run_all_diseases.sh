#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
PIPELINE="$SCRIPT_DIR/run_disease_pipeline.sh"

DISEASES=(
    fabry
    monilethrix
    ais
    nf1
)

MASTER_LOG="$PROJECT_ROOT/results/all_diseases.pipeline.log"

exec > >(tee -a "$MASTER_LOG") 2>&1

echo "============================================================"
echo "STARTING ALL DISEASE PIPELINES"
echo "Time: $(date --iso-8601=seconds)"
echo "============================================================"

for disease in "${DISEASES[@]}"
do
    echo
    echo "############################################################"
    echo "RUNNING: $disease"
    echo "############################################################"

    bash "$PIPELINE" "$disease"

    echo "PASS: $disease completed successfully"
done

echo
echo "============================================================"
echo "ALL FOUR DISEASE PIPELINES COMPLETED SUCCESSFULLY"
echo "============================================================"
echo "Master log: $MASTER_LOG"
