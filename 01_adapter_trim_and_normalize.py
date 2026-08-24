#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Bio-STRATA amplicon preprocessing pipeline
==========================================

This script applies the same preprocessing workflow to multiple R1 CSV files.
By default, it automatically discovers files such as:
    A_R1.csv, B_R1.csv, C_R1.csv, D_R1.csv, F_R1.csv

For each sample, it performs:
1. Similarity-based detection and trimming of the Illumina-compatible adapter
   from the read end.
2. Classification of each 120-bp insert as forward or reverse relative to the
   reference sequence.
3. Reverse-complementation of reverse-oriented inserts so that all retained
   sequences are normalized to the forward orientation.

Per-sample outputs:
    <sample>_R1_step1_adapter_trimmed.xlsx
    <sample>_R1_step2_direction_normalized.xlsx
    <sample>_R1_step3_all_forward.xlsx

Combined output:
    Bio_STRATA_preprocessing_summary.csv

Requirements:
    pip install pandas openpyxl

Examples:
    # Auto-discover all *_R1.csv files in the current folder
    python 01_adapter_trim_and_normalize.py

    # Process selected samples only
    python 01_adapter_trim_and_normalize.py --samples A B C D F

    # Process files in another folder
    python 01_adapter_trim_and_normalize.py --folder /path/to/data
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


# ============================================================
# Analysis settings preserved from the final study workflow
# ============================================================

SEQUENCE_COLUMN = "Sequence"

REFERENCE_FORWARD = (
    "ACGTTGACCTGATCGTACGAGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGC"
    "ACCTACGTTCAATATTACAGGCGAACATACTTACTATGCTAGCTAGGCTACGATCG"
)

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCCAG"

EXPECTED_INSERT_LENGTH = len(REFERENCE_FORWARD)  # 120 bp
MIN_ADAPTER_OVERLAP = 12
MAX_ADAPTER_ERROR_RATE = 0.20
BOUNDARY_SEARCH_MARGIN = 15
MAX_REFERENCE_MISMATCHES = 10


# ============================================================
# Utility functions
# ============================================================

def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def hamming_distance(a: str, b: str) -> int:
    """Return the number of mismatches between equal-length strings."""
    return sum(x != y for x, y in zip(a, b))


def find_best_adapter_boundary(seq: str) -> Dict[str, object]:
    """Find the best adapter-start position near the expected insert boundary."""
    seq = seq.upper().strip()
    seq_len = len(seq)

    start_min = max(
        MIN_ADAPTER_OVERLAP,
        EXPECTED_INSERT_LENGTH - BOUNDARY_SEARCH_MARGIN,
    )
    start_max = min(
        seq_len - MIN_ADAPTER_OVERLAP,
        EXPECTED_INSERT_LENGTH + BOUNDARY_SEARCH_MARGIN,
    )

    candidates: List[Dict[str, object]] = []

    for boundary in range(start_min, start_max + 1):
        observed = seq[boundary:]
        overlap = min(len(observed), len(ADAPTER))

        if overlap < MIN_ADAPTER_OVERLAP:
            continue

        obs_part = observed[:overlap]
        adapter_part = ADAPTER[:overlap]

        mismatches = hamming_distance(obs_part, adapter_part)
        error_rate = mismatches / overlap

        candidates.append(
            {
                "boundary": boundary,
                "overlap": overlap,
                "mismatches": mismatches,
                "error_rate": error_rate,
                "distance_from_expected": abs(boundary - EXPECTED_INSERT_LENGTH),
            }
        )

    if not candidates:
        return {
            "status": "Adapter_not_found",
            "boundary": -1,
            "overlap": 0,
            "mismatches": 0,
            "error_rate": 1.0,
            "trimmed_insert": "",
            "adapter_observed": "",
        }

    best = sorted(
        candidates,
        key=lambda x: (
            x["error_rate"],
            x["mismatches"],
            x["distance_from_expected"],
            -x["overlap"],
        ),
    )[0]

    if float(best["error_rate"]) <= MAX_ADAPTER_ERROR_RATE:
        status = "Adapter_confident"
        boundary = int(best["boundary"])
        trimmed_insert = seq[:boundary]
        adapter_observed = seq[boundary:]
    else:
        status = "Adapter_not_found"
        boundary = -1
        trimmed_insert = ""
        adapter_observed = ""

    return {
        "status": status,
        "boundary": boundary,
        "overlap": best["overlap"],
        "mismatches": best["mismatches"],
        "error_rate": best["error_rate"],
        "trimmed_insert": trimmed_insert,
        "adapter_observed": adapter_observed,
    }


def classify_orientation(insert_seq: str) -> Dict[str, object]:
    """Classify a 120-bp insert as forward, reverse, or unclassified."""
    forward_ref = REFERENCE_FORWARD
    reverse_ref = reverse_complement(REFERENCE_FORWARD)

    if len(insert_seq) != EXPECTED_INSERT_LENGTH:
        return {
            "orientation": "Unclassified",
            "forward_mismatches": None,
            "reverse_mismatches": None,
            "reason": f"Unexpected_insert_length_{len(insert_seq)}",
        }

    f_mm = hamming_distance(insert_seq, forward_ref)
    r_mm = hamming_distance(insert_seq, reverse_ref)

    if min(f_mm, r_mm) > MAX_REFERENCE_MISMATCHES:
        return {
            "orientation": "Unclassified",
            "forward_mismatches": f_mm,
            "reverse_mismatches": r_mm,
            "reason": "Too_many_reference_mismatches",
        }

    if f_mm < r_mm:
        orientation = "Forward"
        reason = "Closer_to_forward_reference"
    elif r_mm < f_mm:
        orientation = "Reverse"
        reason = "Closer_to_reverse_reference"
    else:
        orientation = "Unclassified"
        reason = "Equal_forward_reverse_score"

    return {
        "orientation": orientation,
        "forward_mismatches": f_mm,
        "reverse_mismatches": r_mm,
        "reason": reason,
    }


def autosize_excel_columns(writer) -> None:
    for worksheet in writer.book.worksheets:
        for column_cells in worksheet.columns:
            max_length = 0
            column_letter = column_cells[0].column_letter

            for cell in column_cells[:200]:
                value = "" if cell.value is None else str(cell.value)
                max_length = max(max_length, len(value))

            worksheet.column_dimensions[column_letter].width = min(max_length + 2, 60)

        worksheet.freeze_panes = "A2"


def preferred_input_for_sample(folder: Path, sample: str) -> Optional[Path]:
    """Return the preferred CSV for a sample, supporting legacy '(1)' filenames."""
    exact = folder / f"{sample}_R1.csv"
    if exact.exists():
        return exact

    candidates = sorted(folder.glob(f"{sample}_R1*.csv"))
    return candidates[0] if candidates else None


def discover_samples(folder: Path) -> List[str]:
    """Discover sample labels from files matching <sample>_R1*.csv."""
    samples = set()
    pattern = re.compile(r"^(.+?)_R1(?:\(\d+\))?\.csv$", flags=re.IGNORECASE)

    for path in folder.iterdir():
        if not path.is_file():
            continue
        match = pattern.match(path.name)
        if match:
            samples.add(match.group(1))

    return sorted(samples)


def existing_columns(df: pd.DataFrame, columns: List[str]) -> List[str]:
    """Keep only columns present in a DataFrame, preserving requested order."""
    return [column for column in columns if column in df.columns]


# ============================================================
# Per-sample analysis
# ============================================================

def process_sample(sample: str, input_path: Path, output_folder: Path) -> Dict[str, object]:
    df = pd.read_csv(input_path)

    if SEQUENCE_COLUMN not in df.columns:
        raise ValueError(
            f"{input_path.name}: '{SEQUENCE_COLUMN}' column is missing. "
            f"Available columns: {list(df.columns)}"
        )

    df = df.copy()
    df.insert(0, "Original_row", range(2, len(df) + 2))

    output_records = []

    for _, row in df.iterrows():
        raw_seq = str(row[SEQUENCE_COLUMN]).strip().upper()
        adapter_result = find_best_adapter_boundary(raw_seq)

        orientation_result = {
            "orientation": "Unclassified",
            "forward_mismatches": None,
            "reverse_mismatches": None,
            "reason": "Adapter_not_confident",
        }

        normalized_seq = ""
        normalization_action = "Excluded"

        if adapter_result["status"] == "Adapter_confident":
            insert_seq = str(adapter_result["trimmed_insert"])
            orientation_result = classify_orientation(insert_seq)

            if orientation_result["orientation"] == "Forward":
                normalized_seq = insert_seq
                normalization_action = "Kept_as_forward"
            elif orientation_result["orientation"] == "Reverse":
                normalized_seq = reverse_complement(insert_seq)
                normalization_action = "Reverse_complemented"

        rec = row.to_dict()
        rec.update(
            {
                "Sample": sample,
                "Raw_sequence": raw_seq,
                "Raw_length": len(raw_seq),
                "Adapter_status": adapter_result["status"],
                "Adapter_start_0based": adapter_result["boundary"],
                "Adapter_overlap_checked": adapter_result["overlap"],
                "Adapter_mismatches": adapter_result["mismatches"],
                "Adapter_error_rate": adapter_result["error_rate"],
                "Adapter_observed": adapter_result["adapter_observed"],
                "Trimmed_insert": adapter_result["trimmed_insert"],
                "Trimmed_insert_length": len(str(adapter_result["trimmed_insert"])),
                "Orientation": orientation_result["orientation"],
                "Forward_reference_mismatches": orientation_result["forward_mismatches"],
                "Reverse_reference_mismatches": orientation_result["reverse_mismatches"],
                "Orientation_reason": orientation_result["reason"],
                "Normalization_action": normalization_action,
                "Normalized_forward_sequence": normalized_seq,
                "Normalized_length": len(normalized_seq),
            }
        )
        output_records.append(rec)

    result = pd.DataFrame(output_records)

    output_step1 = output_folder / f"{sample}_R1_step1_adapter_trimmed.xlsx"
    output_step2 = output_folder / f"{sample}_R1_step2_direction_normalized.xlsx"
    output_step3 = output_folder / f"{sample}_R1_step3_all_forward.xlsx"

    # STEP 1
    step1_summary = pd.DataFrame(
        {
            "Metric": [
                "Total input reads",
                "Adapter confidently detected",
                "Adapter not found/confident",
                f"Trimmed insert length = {EXPECTED_INSERT_LENGTH} bp",
                f"Trimmed insert length != {EXPECTED_INSERT_LENGTH} bp",
            ],
            "Value": [
                len(result),
                (result["Adapter_status"] == "Adapter_confident").sum(),
                (result["Adapter_status"] != "Adapter_confident").sum(),
                (result["Trimmed_insert_length"] == EXPECTED_INSERT_LENGTH).sum(),
                (
                    (result["Adapter_status"] == "Adapter_confident")
                    & (result["Trimmed_insert_length"] != EXPECTED_INSERT_LENGTH)
                ).sum(),
            ],
        }
    )

    step1_cols = existing_columns(
        result,
        [
            "Sample",
            "Original_row",
            "ID",
            "Description",
            "Raw_sequence",
            "Raw_length",
            "Adapter_status",
            "Adapter_start_0based",
            "Adapter_overlap_checked",
            "Adapter_mismatches",
            "Adapter_error_rate",
            "Adapter_observed",
            "Trimmed_insert",
            "Trimmed_insert_length",
        ],
    )

    with pd.ExcelWriter(output_step1, engine="openpyxl") as writer:
        step1_summary.to_excel(writer, sheet_name="Summary", index=False)
        result[step1_cols].to_excel(writer, sheet_name="Adapter_trimmed", index=False)
        result[result["Adapter_status"] != "Adapter_confident"].to_excel(
            writer, sheet_name="Adapter_failed", index=False
        )
        autosize_excel_columns(writer)

    # STEP 2
    direction_counts = result.groupby("Orientation").size().reset_index(name="Read_count")
    direction_counts["Percent"] = direction_counts["Read_count"] / len(result) * 100

    step2_cols = existing_columns(
        result,
        [
            "Sample",
            "Original_row",
            "ID",
            "Description",
            "Trimmed_insert",
            "Trimmed_insert_length",
            "Orientation",
            "Forward_reference_mismatches",
            "Reverse_reference_mismatches",
            "Orientation_reason",
            "Normalization_action",
            "Normalized_forward_sequence",
            "Normalized_length",
        ],
    )

    with pd.ExcelWriter(output_step2, engine="openpyxl") as writer:
        direction_counts.to_excel(writer, sheet_name="Summary", index=False)
        result[step2_cols].to_excel(writer, sheet_name="Direction_normalized", index=False)
        result[result["Orientation"] == "Unclassified"].to_excel(
            writer, sheet_name="Unclassified", index=False
        )
        autosize_excel_columns(writer)

    # STEP 3
    included = result[
        result["Orientation"].isin(["Forward", "Reverse"])
        & result["Normalized_forward_sequence"].ne("")
    ].copy()

    excluded = result.loc[~result.index.isin(included.index)].copy()

    forward_cols = existing_columns(
        included,
        [
            "Sample",
            "Original_row",
            "ID",
            "Description",
            "Orientation",
            "Normalization_action",
            "Normalized_forward_sequence",
        ],
    )

    all_forward_reads = included[forward_cols].rename(
        columns={"Normalized_forward_sequence": "Sequence"}
    )

    unique_sequences = (
        all_forward_reads.groupby("Sequence")
        .size()
        .reset_index(name="Count")
        .sort_values(["Count", "Sequence"], ascending=[False, True])
        .reset_index(drop=True)
    )
    unique_sequences.insert(0, "Rank", range(1, len(unique_sequences) + 1))

    step3_summary = pd.DataFrame(
        {
            "Metric": [
                "Total input reads",
                "Included all-forward reads",
                "Original forward reads",
                "Reverse reads converted",
                "Excluded reads",
                "Unique normalized sequences",
                f"All normalized sequence length = {EXPECTED_INSERT_LENGTH} bp",
            ],
            "Value": [
                len(result),
                len(included),
                (included["Orientation"] == "Forward").sum(),
                (included["Orientation"] == "Reverse").sum(),
                len(excluded),
                len(unique_sequences),
                (included["Normalized_length"] == EXPECTED_INSERT_LENGTH).all(),
            ],
        }
    )

    with pd.ExcelWriter(output_step3, engine="openpyxl") as writer:
        step3_summary.to_excel(writer, sheet_name="Summary", index=False)
        all_forward_reads.to_excel(writer, sheet_name="All_forward_reads", index=False)
        unique_sequences.to_excel(writer, sheet_name="Grouped_sequences", index=False)
        included.to_excel(writer, sheet_name="All_forward_details", index=False)
        excluded.to_excel(writer, sheet_name="Excluded_reads", index=False)
        autosize_excel_columns(writer)

    return {
        "sample": sample,
        "input_file": input_path.name,
        "total_input_reads": len(result),
        "adapter_confident_reads": int((result["Adapter_status"] == "Adapter_confident").sum()),
        "included_all_forward_reads": len(included),
        "original_forward_reads": int((included["Orientation"] == "Forward").sum()),
        "reverse_reads_converted": int((included["Orientation"] == "Reverse").sum()),
        "excluded_reads": len(excluded),
        "unique_normalized_sequences": len(unique_sequences),
        "output_step3": output_step3.name,
    }


# ============================================================
# Main
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Preprocess Bio-STRATA R1 amplicon CSV files for one or more samples."
    )
    parser.add_argument(
        "--folder",
        type=Path,
        default=Path.cwd(),
        help="Folder containing *_R1.csv files. Default: current folder.",
    )
    parser.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Optional sample labels to process, e.g. --samples A B C D F. "
        "If omitted, samples are auto-discovered.",
    )
    args = parser.parse_args()

    folder = args.folder.resolve()
    if not folder.exists() or not folder.is_dir():
        raise SystemExit(f"Folder not found: {folder}")

    samples = args.samples if args.samples else discover_samples(folder)
    if not samples:
        raise SystemExit("No *_R1.csv input files were found.")

    print(f"Folder: {folder}")
    print(f"Samples: {', '.join(samples)}")
    print(f"Expected insert length: {EXPECTED_INSERT_LENGTH} bp")
    print()

    summaries = []

    for sample in samples:
        input_path = preferred_input_for_sample(folder, sample)
        if input_path is None:
            print(f"[SKIP] {sample}: no matching {sample}_R1*.csv file found")
            continue

        print(f"Processing {sample}: {input_path.name}")
        summary = process_sample(sample, input_path, folder)
        summaries.append(summary)
        print(
            f"  included={summary['included_all_forward_reads']:,} / "
            f"{summary['total_input_reads']:,}; "
            f"excluded={summary['excluded_reads']:,}"
        )

    if not summaries:
        raise SystemExit("No samples were processed.")

    summary_path = folder / "Bio_STRATA_preprocessing_summary.csv"
    pd.DataFrame(summaries).to_csv(summary_path, index=False, encoding="utf-8-sig")

    print()
    print("Completed.")
    print(f"Combined summary: {summary_path.name}")


if __name__ == "__main__":
    main()
