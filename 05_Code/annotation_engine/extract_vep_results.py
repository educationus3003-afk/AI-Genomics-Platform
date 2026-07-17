#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


def extract_vep(vcf_path: Path, output_path: Path) -> None:
    csq_fields: list[str] = []
    rows: list[dict[str, str]] = []

    with vcf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r"Format: ([^\">]+)", line)
                if match:
                    csq_fields = match.group(1).split("|")

            elif line.startswith("#"):
                continue

            else:
                columns = line.rstrip("\n").split("\t")

                if len(columns) < 8:
                    continue

                chrom, pos, variant_id, ref, alt, _, _, info = columns[:8]

                info_values = {}
                for item in info.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        info_values[key] = value

                csq_entries = info_values.get("CSQ", "").split(",")

                for entry in csq_entries:
                    if not entry:
                        continue

                    values = entry.split("|")
                    annotation = dict(zip(csq_fields, values))

                    rows.append(
                        {
                            "CHROM": chrom,
                            "POS": pos,
                            "ID": variant_id,
                            "REF": ref,
                            "ALT": alt,
                            "Gene": annotation.get("SYMBOL", ""),
                            "Consequence": annotation.get("Consequence", ""),
                            "Impact": annotation.get("IMPACT", ""),
                            "Feature": annotation.get("Feature", ""),
                            "HGVSc": annotation.get("HGVSc", ""),
                            "HGVSp": annotation.get("HGVSp", ""),
                            "Existing_variation": annotation.get(
                                "Existing_variation", ""
                            ),
                            "Canonical": annotation.get("CANONICAL", ""),
                            "MANE_Select": annotation.get("MANE_SELECT", ""),
                            "ClinVar_significance": annotation.get(
                                "CLIN_SIG", ""
                            ),
                            "gnomAD_exome_AF": annotation.get(
                                "gnomADe_AF", ""
                            ),
                            "gnomAD_genome_AF": annotation.get(
                                "gnomADg_AF", ""
                            ),
                        }
                    )

    fieldnames = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "Gene",
        "Consequence",
        "Impact",
        "Feature",
        "HGVSc",
        "HGVSp",
        "Existing_variation",
        "Canonical",
        "MANE_Select",
        "ClinVar_significance",
        "gnomAD_exome_AF",
        "gnomAD_genome_AF",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(
            "Usage: extract_vep_results.py input.vep.vcf output.tsv"
        )

    extract_vep(Path(sys.argv[1]), Path(sys.argv[2]))


if __name__ == "__main__":
    main()
