#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

MASTER="$ROOT/08_Reports/AI_Genomics_Platform_Complete_Report.html"
PAGES="$ROOT/docs/index.html"

if [[ ! -f "$MASTER" ]]; then
    echo "ERROR: Master report does not exist:"
    echo "$MASTER"
    exit 1
fi

mkdir -p "$ROOT/docs"

cp "$MASTER" "$PAGES"

echo
echo "Report synchronized successfully."
echo
sha256sum "$MASTER" "$PAGES"
