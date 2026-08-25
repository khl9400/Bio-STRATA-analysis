#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Bio-STRATA R2 SNP and MT-fraction analysis
==========================================

This script applies the same final SNP-calling workflow used for R1
to all normalized R2 outputs.

Input:
    <sample>_R2_step3_all_forward.xlsx

Output:
    <sample>_R2_SNP_result.xlsx

Combined output:
    Bio_STRATA_R2_SNP_summary.csv

Default samples:
    A, B, C, D, F

SNP definition:
    Position: 52 (1-based)
    WT allele: T
    MT allele: C
    Other base / insufficient length: Unknown

MT fraction:
    MT / (WT + MT) * 100

Requirements:
    pip install pandas openpyxl
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd


# ============================================================
# Analysis settings
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

            worksheet.column_dimensions[column_letter].width = min(
                max_length + 2,
                65
            )


# ============================================================
# Per-sample analysis
# ============================================================

def analyze_sample(
    sample: str,
    input_path: Path,
    output_folder: Path,
) -> Dict[str, object]:

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

    if result.empty:
        raise ValueError(
            f"{input_path.name}: no reads found in '{INPUT_SHEET}'."
        )

    wt_count = int((result["Classification"] == "WT").sum())
    mt_count = int((result["Classification"] == "MT").sum())
    unknown_count = int((result["Classification"] == "Unknown").sum())
    valid_count = wt_count + mt_count

    wt_ratio = (
        wt_count / valid_count * 100
        if valid_count > 0
        else float("nan")
    )
    mt_ratio = (
        mt_count / valid_count * 100
        if valid_count > 0
        else float("nan")
    )

    # EXACTLY the same Summary layout as the R1 result.
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

    output_file = output_folder / f"{sample}_R2_SNP_result.xlsx"

    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        classification_counts.to_excel(
            writer,
            sheet_name="Classification_counts",
            index=False,
        )
        base_counts.to_excel(
            writer,
            sheet_name="SNP_base_counts",
            index=False,
        )
        result.to_excel(
            writer,
            sheet_name="Read_results",
            index=False,
        )
        result[result["Classification"] == "MT"].to_excel(
            writer,
            sheet_name="MT_reads",
            index=False,
        )
        result[result["Classification"] == "Unknown"].to_excel(
            writer,
            sheet_name="Unknown_reads",
            index=False,
        )
        autosize_excel_columns(writer)

    # Print the same final result table to the terminal.
    print()
    print(f"Sample: {sample}")
    print(f"Input file: {input_path.name}")
    print(f"Output file: {output_file.name}")
    print("=" * 60)
    print(summary.to_string(index=False))

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
        description=(
            "Calculate WT/MT/Unknown read counts and MT fractions "
            "for Bio-STRATA R2 samples."
        )
    )
    parser.add_argument(
        "--folder",
        type=Path,
        default=Path.cwd(),
        help="Folder containing *_R2_step3_all_forward.xlsx files.",
    )
    parser.add_argument(
        "--samples",
        nargs="*",
        default=["A", "B", "C", "D", "F"],
        help="Samples to analyze. Default: A B C D F",
    )
    args = parser.parse_args()

    folder = args.folder.resolve()

    if not folder.exists() or not folder.is_dir():
        raise SystemExit(f"Folder not found: {folder}")

    summaries: List[Dict[str, object]] = []

    for sample in args.samples:
        input_path = folder / f"{sample}_R2_step3_all_forward.xlsx"

        if not input_path.exists():
            print(
                f"[SKIP] {sample}: "
                f"{input_path.name} not found"
            )
            continue

        summary = analyze_sample(
            sample=sample,
            input_path=input_path,
            output_folder=folder,
        )
        summaries.append(summary)

    if not summaries:
        raise SystemExit("No R2 samples were analyzed.")

    combined_path = folder / "Bio_STRATA_R2_SNP_summary.csv"
    pd.DataFrame(summaries).to_csv(
        combined_path,
        index=False,
        encoding="utf-8-sig",
    )

    print()
    print("Completed.")
    print(f"Combined summary: {combined_path.name}")


if __name__ == "__main__":
    main()
