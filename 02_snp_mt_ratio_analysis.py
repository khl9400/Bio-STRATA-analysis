#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Bio-STRATA SNP and MT-fraction analysis
=======================================

This script applies the final SNP-calling rule to all normalized R1 outputs
created by 01_adapter_trim_and_normalize.py.

Default SNP definition preserved from the final study workflow:
    SNP position: 52 (1-based)
    WT allele: T
    MT allele: C
    Other base or insufficient sequence length: Unknown

MT fraction:
    MT / (WT + MT) * 100

Input files are automatically discovered as:
    <sample>_R1_step3_all_forward.xlsx

Per-sample output:
    <sample>_R1_SNP_result.xlsx

Combined output:
    Bio_STRATA_SNP_summary.csv

Requirements:
    pip install pandas openpyxl

Examples:
    # Auto-discover all normalized samples in the current folder
    python 02_snp_mt_ratio_analysis.py

    # Analyze selected samples only
    python 02_snp_mt_ratio_analysis.py --samples A B C D F

    # Analyze files in another folder
    python 02_snp_mt_ratio_analysis.py --folder /path/to/data
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List

import pandas as pd


# ============================================================
# Analysis settings preserved from the final study workflow
# ============================================================

INPUT_SHEET = "All_forward_reads"
SEQUENCE_COLUMN = "Sequence"

SNP_POSITION = 52  # 1-based
WT_BASE = "T"
MT_BASE = "C"

REFERENCE = (
    "ACGTTGACCTGATCGTACGAGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGC"
    "ACCTACGTTCAATATTACAGGCGAACATACTTACTATGCTAGCTAGGCTACGATCG"
)


# ============================================================
# Utility functions
# ============================================================

def count_reference_mismatches(sequence: str) -> int | None:
    """Count mismatches against the full reference when lengths match."""
    if len(sequence) != len(REFERENCE):
        return None
    return sum(a != b for a, b in zip(sequence, REFERENCE))


def autosize_excel_columns(writer) -> None:
    for worksheet in writer.book.worksheets:
        worksheet.freeze_panes = "A2"

        for column_cells in worksheet.columns:
            column_letter = column_cells[0].column_letter
            max_length = 0

            for cell in column_cells[:300]:
                value = "" if cell.value is None else str(cell.value)
                max_length = max(max_length, len(value))

            worksheet.column_dimensions[column_letter].width = min(max_length + 2, 65)


def discover_inputs(folder: Path) -> Dict[str, Path]:
    """Discover <sample>_R1_step3_all_forward.xlsx files."""
    found: Dict[str, Path] = {}
    pattern = re.compile(r"^(.+?)_R1_step3_all_forward\.xlsx$", flags=re.IGNORECASE)

    for path in folder.iterdir():
        if not path.is_file():
            continue
        match = pattern.match(path.name)
        if match:
            found[match.group(1)] = path

    return dict(sorted(found.items()))


# ============================================================
# Per-sample SNP analysis
# ============================================================

def analyze_sample(sample: str, input_path: Path, output_folder: Path) -> Dict[str, object]:
    df = pd.read_excel(input_path, sheet_name=INPUT_SHEET)

    if SEQUENCE_COLUMN not in df.columns:
        raise ValueError(
            f"{input_path.name}: '{SEQUENCE_COLUMN}' column is missing. "
            f"Available columns: {list(df.columns)}"
        )

    snp_index = SNP_POSITION - 1
    records = []

    for _, row in df.iterrows():
        sequence = str(row[SEQUENCE_COLUMN]).strip().upper()

        if len(sequence) <= snp_index:
            snp_base = ""
            classification = "Unknown"
            reason = "Sequence_too_short"
        else:
            snp_base = sequence[snp_index]

            if snp_base == WT_BASE:
                classification = "WT"
                reason = "WT_base_at_SNP"
            elif snp_base == MT_BASE:
                classification = "MT"
                reason = "MT_base_at_SNP"
            else:
                classification = "Unknown"
                reason = f"Unexpected_base_{snp_base}"

        record = row.to_dict()
        record.update(
            {
                "Sample": sample,
                "Sequence_length": len(sequence),
                "Reference_mismatches": count_reference_mismatches(sequence),
                "SNP_position_1based": SNP_POSITION,
                "SNP_base": snp_base,
                "Classification": classification,
                "Classification_reason": reason,
            }
        )
        records.append(record)

    result = pd.DataFrame(records)

    wt_count = int((result["Classification"] == "WT").sum())
    mt_count = int((result["Classification"] == "MT").sum())
    unknown_count = int((result["Classification"] == "Unknown").sum())
    valid_count = wt_count + mt_count

    mt_ratio = mt_count / valid_count * 100 if valid_count > 0 else float("nan")
    wt_ratio = wt_count / valid_count * 100 if valid_count > 0 else float("nan")

    summary = pd.DataFrame(
        {
            "Metric": [
                "Sample",
                "Total reads",
                "Valid allele reads (WT + MT)",
                "WT reads",
                "MT reads",
                "Unknown reads",
                "WT ratio (%)",
                "MT ratio (%)",
                "SNP position (1-based)",
                "WT base",
                "MT base",
            ],
            "Value": [
                sample,
                len(result),
                valid_count,
                wt_count,
                mt_count,
                unknown_count,
                round(wt_ratio, 6),
                round(mt_ratio, 6),
                SNP_POSITION,
                WT_BASE,
                MT_BASE,
            ],
        }
    )

    base_counts = (
        result["SNP_base"]
        .replace("", "No_base")
        .value_counts(dropna=False)
        .rename_axis("SNP_base")
        .reset_index(name="Read_count")
    )
    base_counts["Percent_of_total"] = (
        base_counts["Read_count"] / len(result) * 100
    ).round(6)

    classification_counts = (
        result["Classification"]
        .value_counts()
        .rename_axis("Classification")
        .reset_index(name="Read_count")
    )
    classification_counts["Percent_of_total"] = (
        classification_counts["Read_count"] / len(result) * 100
    ).round(6)

    output_file = output_folder / f"{sample}_R1_SNP_result.xlsx"

    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        classification_counts.to_excel(writer, sheet_name="Classification_counts", index=False)
        base_counts.to_excel(writer, sheet_name="SNP_base_counts", index=False)
        result.to_excel(writer, sheet_name="Read_results", index=False)
        result[result["Classification"] == "MT"].to_excel(
            writer, sheet_name="MT_reads", index=False
        )
        result[result["Classification"] == "Unknown"].to_excel(
            writer, sheet_name="Unknown_reads", index=False
        )
        autosize_excel_columns(writer)

    return {
        "sample": sample,
        "input_file": input_path.name,
        "total_reads": len(result),
        "valid_allele_reads": valid_count,
        "WT_reads": wt_count,
        "MT_reads": mt_count,
        "Unknown_reads": unknown_count,
        "WT_ratio_pct": round(wt_ratio, 6),
        "MT_ratio_pct": round(mt_ratio, 6),
        "SNP_position_1based": SNP_POSITION,
        "WT_base": WT_BASE,
        "MT_base": MT_BASE,
        "output_file": output_file.name,
    }


# ============================================================
# Main
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate WT/MT/Unknown read counts and MT fractions for Bio-STRATA samples."
    )
    parser.add_argument(
        "--folder",
        type=Path,
        default=Path.cwd(),
        help="Folder containing *_R1_step3_all_forward.xlsx files. Default: current folder.",
    )
    parser.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Optional sample labels to analyze, e.g. --samples A B C D F. "
        "If omitted, samples are auto-discovered.",
    )
    args = parser.parse_args()

    folder = args.folder.resolve()
    if not folder.exists() or not folder.is_dir():
        raise SystemExit(f"Folder not found: {folder}")

    discovered = discover_inputs(folder)
    if not discovered:
        raise SystemExit("No *_R1_step3_all_forward.xlsx files were found.")

    if args.samples:
        selected = {sample: discovered[sample] for sample in args.samples if sample in discovered}
        missing = [sample for sample in args.samples if sample not in discovered]
        for sample in missing:
            print(f"[SKIP] {sample}: {sample}_R1_step3_all_forward.xlsx not found")
    else:
        selected = discovered

    if not selected:
        raise SystemExit("No requested samples were available for analysis.")

    print(f"Folder: {folder}")
    print(f"Samples: {', '.join(selected)}")
    print(f"SNP position: {SNP_POSITION} (1-based)")
    print(f"WT/MT: {WT_BASE}/{MT_BASE}")
    print()

    summaries: List[Dict[str, object]] = []

    for sample, input_path in selected.items():
        print(f"Analyzing {sample}: {input_path.name}")
        summary = analyze_sample(sample, input_path, folder)
        summaries.append(summary)
        print(
            f"  WT={summary['WT_reads']:,}; MT={summary['MT_reads']:,}; "
            f"Unknown={summary['Unknown_reads']:,}; "
            f"MT fraction={summary['MT_ratio_pct']:.6f}%"
        )

    combined_summary = pd.DataFrame(summaries)
    combined_path = folder / "Bio_STRATA_SNP_summary.csv"
    combined_summary.to_csv(combined_path, index=False, encoding="utf-8-sig")

    print()
    print("Completed.")
    print(f"Combined summary: {combined_path.name}")


if __name__ == "__main__":
    main()
