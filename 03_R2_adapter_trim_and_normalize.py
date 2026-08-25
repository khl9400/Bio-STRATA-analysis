#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A~F R2.csv 전처리 파이프라인

기존 A_R2 전처리 코드의 분석 조건과 처리 로직은 그대로 유지하고,
A_R2.csv, B_R2.csv, C_R2.csv, D_R2.csv, F_R2.csv에 동일하게 적용하도록 확장함.

각 sample에 대해:
1) read 말단에서 Illumina adapter를 유사도 기반으로 탐지하여 제거
2) adapter가 제거된 insert가 Forward인지 Reverse인지 판별
3) Reverse insert는 reverse-complement하여 모두 Forward 방향으로 통일

출력 예:
- A_R2_step1_adapter_trimmed.xlsx
- A_R2_step2_direction_normalized.xlsx
- A_R2_step3_all_forward.xlsx
- B_R2_step1_adapter_trimmed.xlsx
...

필요 패키지:
    pip install pandas openpyxl
"""

from pathlib import Path
import argparse
import pandas as pd


# ============================================================
# 사용자 설정
# ============================================================

SEQUENCE_COLUMN = "Sequence"

REFERENCE_FORWARD = (
    "ACGTTGACCTGATCGTACGAGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGC"
    "ACCTACGTTCAATATTACAGGCGAACATACTTACTATGCTAGCTAGGCTACGATCG"
)

ADAPTER = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

EXPECTED_INSERT_LENGTH = len(REFERENCE_FORWARD)  # 120 bp

# adapter 탐지 시 최소로 확인할 adapter overlap
MIN_ADAPTER_OVERLAP = 12

# adapter overlap 안에서 허용할 최대 mismatch 비율
MAX_ADAPTER_ERROR_RATE = 0.20

# 예상 insert 길이 근처에서 adapter를 우선 탐색할 범위
BOUNDARY_SEARCH_MARGIN = 15

# 방향 판별 시 reference 대비 허용할 최대 mismatch 수
MAX_REFERENCE_MISMATCHES = 10


# ============================================================
# 기본 함수
# ============================================================

def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def hamming_distance(a: str, b: str) -> int:
    """길이가 같은 문자열의 mismatch 수."""
    return sum(x != y for x, y in zip(a, b))


def mismatch_rate(a: str, b: str) -> float:
    if len(a) == 0:
        return 1.0
    return hamming_distance(a, b) / len(a)


def find_best_adapter_boundary(seq: str):
    """
    read의 뒤쪽에서 adapter prefix와 가장 잘 맞는 시작점을 찾음.

    반환:
        {
            status,
            boundary,
            overlap,
            mismatches,
            error_rate,
            trimmed_insert,
            adapter_observed
        }
    """

    seq = seq.upper().strip()
    seq_len = len(seq)

    # adapter는 insert 뒤쪽에 있어야 하므로 예상 경계 주변을 우선 탐색
    start_min = max(
        MIN_ADAPTER_OVERLAP,
        EXPECTED_INSERT_LENGTH - BOUNDARY_SEARCH_MARGIN
    )
    start_max = min(
        seq_len - MIN_ADAPTER_OVERLAP,
        EXPECTED_INSERT_LENGTH + BOUNDARY_SEARCH_MARGIN
    )

    candidates = []

    for boundary in range(start_min, start_max + 1):
        observed = seq[boundary:]
        overlap = min(len(observed), len(ADAPTER))

        if overlap < MIN_ADAPTER_OVERLAP:
            continue

        obs_part = observed[:overlap]
        adapter_part = ADAPTER[:overlap]

        mismatches = hamming_distance(obs_part, adapter_part)
        error_rate = mismatches / overlap

        candidates.append({
            "boundary": boundary,
            "overlap": overlap,
            "mismatches": mismatches,
            "error_rate": error_rate,
            "distance_from_expected": abs(boundary - EXPECTED_INSERT_LENGTH),
        })

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

    # 1) error rate 최소
    # 2) mismatch 수 최소
    # 3) 예상 insert 길이 120 bp에 가까운 경계 우선
    best = sorted(
        candidates,
        key=lambda x: (
            x["error_rate"],
            x["mismatches"],
            x["distance_from_expected"],
            -x["overlap"],
        )
    )[0]

    if best["error_rate"] <= MAX_ADAPTER_ERROR_RATE:
        status = "Adapter_confident"
        boundary = best["boundary"]
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


def classify_orientation(insert_seq: str):
    """
    insert 전체를 Forward reference와 Reverse-complement reference에 비교하여
    더 가까운 방향을 선택.
    """

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


def autosize_excel_columns(writer):
    """각 시트의 열 너비를 간단히 조정."""
    for worksheet in writer.book.worksheets:
        for column_cells in worksheet.columns:
            max_length = 0
            column_letter = column_cells[0].column_letter

            for cell in column_cells[:200]:
                value = "" if cell.value is None else str(cell.value)
                max_length = max(max_length, len(value))

            worksheet.column_dimensions[column_letter].width = min(max_length + 2, 60)
        worksheet.freeze_panes = "A2"


# ============================================================
# sample별 실행
# ============================================================

def find_input_file(sample: str):
    candidates = [
        f"{sample}_R2.csv",
        f"{sample}_R2(1).csv",
    ]

    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            return path

    return None


def process_sample(sample: str):

    input_path = find_input_file(sample)

    if input_path is None:
        print(f"[SKIP] {sample}: {sample}_R2.csv 또는 {sample}_R2(1).csv를 찾을 수 없습니다.")
        return

    output_step1 = f"{sample}_R2_step1_adapter_trimmed.xlsx"
    output_step2 = f"{sample}_R2_step2_direction_normalized.xlsx"
    output_step3 = f"{sample}_R2_step3_all_forward.xlsx"

    df = pd.read_csv(input_path)

    if SEQUENCE_COLUMN not in df.columns:
        raise ValueError(
            f"'{SEQUENCE_COLUMN}' 열이 없습니다. 현재 열: {list(df.columns)}"
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

            insert_seq = adapter_result["trimmed_insert"]
            orientation_result = classify_orientation(insert_seq)

            if orientation_result["orientation"] == "Forward":
                normalized_seq = insert_seq
                normalization_action = "Kept_as_forward"

            elif orientation_result["orientation"] == "Reverse":
                normalized_seq = reverse_complement(insert_seq)
                normalization_action = "Reverse_complemented"

        rec = row.to_dict()

        rec.update({
            "Raw_sequence": raw_seq,
            "Raw_length": len(raw_seq),

            "Adapter_status": adapter_result["status"],
            "Adapter_start_0based": adapter_result["boundary"],
            "Adapter_overlap_checked": adapter_result["overlap"],
            "Adapter_mismatches": adapter_result["mismatches"],
            "Adapter_error_rate": adapter_result["error_rate"],
            "Adapter_observed": adapter_result["adapter_observed"],

            "Trimmed_insert": adapter_result["trimmed_insert"],
            "Trimmed_insert_length": len(adapter_result["trimmed_insert"]),

            "Orientation": orientation_result["orientation"],
            "Forward_reference_mismatches":
                orientation_result["forward_mismatches"],
            "Reverse_reference_mismatches":
                orientation_result["reverse_mismatches"],
            "Orientation_reason": orientation_result["reason"],

            "Normalization_action": normalization_action,
            "Normalized_forward_sequence": normalized_seq,
            "Normalized_length": len(normalized_seq),
        })

        output_records.append(rec)

    result = pd.DataFrame(output_records)


    # ========================================================
    # STEP 1: adapter 제거 결과
    # ========================================================

    step1_summary = pd.DataFrame({
        "Metric": [
            "Total input reads",
            "Adapter confidently detected",
            "Adapter not found/confident",
            "Trimmed insert length = 120 bp",
            "Trimmed insert length != 120 bp",
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
    })

    with pd.ExcelWriter(output_step1, engine="openpyxl") as writer:
        step1_summary.to_excel(writer, sheet_name="Summary", index=False)

        result[
            [
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
            ]
        ].to_excel(writer, sheet_name="Adapter_trimmed", index=False)

        result[
            result["Adapter_status"] != "Adapter_confident"
        ].to_excel(writer, sheet_name="Adapter_failed", index=False)

        autosize_excel_columns(writer)


    # ========================================================
    # STEP 2: 방향 판별 및 변환
    # ========================================================

    direction_counts = (
        result.groupby("Orientation")
        .size()
        .reset_index(name="Read_count")
    )

    direction_counts["Percent"] = (
        direction_counts["Read_count"] / len(result) * 100
    )

    with pd.ExcelWriter(output_step2, engine="openpyxl") as writer:
        direction_counts.to_excel(writer, sheet_name="Summary", index=False)

        result[
            [
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
            ]
        ].to_excel(writer, sheet_name="Direction_normalized", index=False)

        result[
            result["Orientation"] == "Unclassified"
        ].to_excel(writer, sheet_name="Unclassified", index=False)

        autosize_excel_columns(writer)


    # ========================================================
    # STEP 3: Forward 방향 통합 파일
    # ========================================================

    included = result[
        result["Orientation"].isin(["Forward", "Reverse"])
        & result["Normalized_forward_sequence"].ne("")
    ].copy()

    excluded = result.loc[~result.index.isin(included.index)].copy()

    all_forward_reads = included[
        [
            "Original_row",
            "ID",
            "Description",
            "Orientation",
            "Normalization_action",
            "Normalized_forward_sequence",
        ]
    ].rename(columns={
        "Normalized_forward_sequence": "Sequence"
    })

    unique_sequences = (
        all_forward_reads.groupby("Sequence")
        .size()
        .reset_index(name="Count")
        .sort_values(
            ["Count", "Sequence"],
            ascending=[False, True]
        )
        .reset_index(drop=True)
    )

    unique_sequences.insert(
        0,
        "Rank",
        range(1, len(unique_sequences) + 1)
    )

    step3_summary = pd.DataFrame({
        "Metric": [
            "Total input reads",
            "Included all-forward reads",
            "Original forward reads",
            "Reverse reads converted",
            "Excluded reads",
            "Unique normalized sequences",
            "All normalized sequence length = 120 bp",
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
    })

    with pd.ExcelWriter(output_step3, engine="openpyxl") as writer:
        step3_summary.to_excel(writer, sheet_name="Summary", index=False)
        all_forward_reads.to_excel(writer, sheet_name="All_forward_reads", index=False)
        unique_sequences.to_excel(writer, sheet_name="Grouped_sequences", index=False)
        included.to_excel(writer, sheet_name="All_forward_details", index=False)
        excluded.to_excel(writer, sheet_name="Excluded_reads", index=False)

        autosize_excel_columns(writer)


    print()
    print("완료")
    print("=" * 60)
    print(f"샘플: {sample}")
    print(f"입력 파일: {input_path}")
    print(f"예상 insert 길이: {EXPECTED_INSERT_LENGTH} bp")
    print(f"1단계: {output_step1}")
    print(f"2단계: {output_step2}")
    print(f"3단계: {output_step3}")
    print("=" * 60)
    print(step3_summary.to_string(index=False))


# ============================================================
# 메인 실행
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--samples",
        nargs="*",
        default=["A", "B", "C", "D", "F"],
        help="분석할 sample 이름. 기본값: A B C D F",
    )
    args = parser.parse_args()

    for sample in args.samples:
        process_sample(sample)


if __name__ == "__main__":
    main()
