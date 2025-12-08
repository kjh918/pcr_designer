# primer/qc.py

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, List, Iterable, Any

import primer3


def compute_heterodimer(f_seq: str, r_seq: str) -> Tuple[float, float]:
    """
    두 올리고 간 heterodimer ΔG / Tm 계산.
    """
    hetero = primer3.calc_heterodimer(f_seq, r_seq)
    if hetero.structure_found:
        het_dg = hetero.dg / 1000.0  # primer3는 1000배 단위
        het_tm = hetero.tm
    else:
        het_dg = 0.0
        het_tm = 0.0
    return het_dg, het_tm


@dataclass
class QCThresholds:
    hairpin_max_tm: float = 47.0
    hairpin_min_dg: float = -5.0
    homodimer_min_dg: float = -6.0
    heterodimer_min_dg: float = -6.0
    heterodimer_max_tm: float = 45.0


def _qc_bool_flags(amp: Dict[str, Any], th: QCThresholds) -> Tuple[bool, bool, bool]:
    """
    내부용: hairpin / homodimer / (FR) heterodimer 각각의 True/False 리턴.
    heterodimer는 기본적으로 F–R 쌍(heterodimer_dg / heterodimer_tm)을 기준으로 판단.
    """
    # Hairpin (Tm, dG 기준)
    hairpin_ok = (
        amp.get("forward_hairpin_tm", 0.0) <= th.hairpin_max_tm
        and amp.get("reverse_hairpin_tm", 0.0) <= th.hairpin_max_tm
        and amp.get("forward_hairpin_dg", 0.0) >= th.hairpin_min_dg
        and amp.get("reverse_hairpin_dg", 0.0) >= th.hairpin_min_dg
    )

    # Homodimer (각각 dG 기준)
    homodimer_ok = (
        amp.get("forward_homodimer_dg", 0.0) >= th.homodimer_min_dg
        and amp.get("reverse_homodimer_dg", 0.0) >= th.homodimer_min_dg
    )

    # F–R heterodimer (기본 heterodimer_dg / heterodimer_tm 사용)
    hetero_fr_ok = (
        amp.get("heterodimer_dg", 0.0) >= th.heterodimer_min_dg
        and amp.get("heterodimer_tm", 0.0) <= th.heterodimer_max_tm
    )

    return hairpin_ok, homodimer_ok, hetero_fr_ok


def amplicon_passes_qc(amp: Dict[str, Any], th: QCThresholds) -> bool:
    """
    최종 QC 통과 여부 (기본적으로 hairpin/homodimer/F–R heterodimer 기준).
    FP/RP까지 포함한 최종 판정은 evaluate_amplicons 안에서 처리.
    """
    hairpin_ok, homodimer_ok, hetero_fr_ok = _qc_bool_flags(amp, th)
    return hairpin_ok and homodimer_ok and hetero_fr_ok


def _hetero_ok(dg: float, tm: float, th: QCThresholds) -> bool:
    """
    단일 heterodimer 쌍에 대한 QC 여부.
    """
    return (dg >= th.heterodimer_min_dg) and (tm <= th.heterodimer_max_tm)


def evaluate_amplicons(
    amplicons: Iterable[Any],
    qc_thresholds: QCThresholds,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Amplicon 객체 리스트에 대해:
      - Amplicon.to_dict() 호출
      - F–R / F–P / R–P heterodimer_dg / heterodimer_tm 계산 추가
      - 각 QC 항목별 결과(O/X) 컬럼 및 최종 QC_PASS 컬럼 추가
      - QC 통과 여부 필터링

    Returns
    -------
    total_rows, filtered_rows : (List[Dict], List[Dict])
        - total_rows: 모든 amplicon의 dict
        - filtered_rows: QC 통과한 amplicon의 dict
    """
    total_rows: List[Dict[str, Any]] = []
    filtered_rows: List[Dict[str, Any]] = []

    for amplicon in amplicons:
        a_dict = amplicon.to_dict()

        f_seq = a_dict.get("forward_sequence")
        r_seq = a_dict.get("reverse_sequence")
        p_seq = a_dict.get("probe_sequence")

        # ---------- 1) heterodimer 계산 ----------
        # F–R
        if f_seq and r_seq:
            het_fr_dg, het_fr_tm = compute_heterodimer(f_seq, r_seq)
        else:
            het_fr_dg, het_fr_tm = 0.0, 0.0

        # F–P
        if f_seq and p_seq:
            het_fp_dg, het_fp_tm = compute_heterodimer(f_seq, p_seq)
        else:
            het_fp_dg, het_fp_tm = 0.0, 0.0

        # R–P
        if r_seq and p_seq:
            het_rp_dg, het_rp_tm = compute_heterodimer(r_seq, p_seq)
        else:
            het_rp_dg, het_rp_tm = 0.0, 0.0

        # 기존 heterodimer_dg / tm 은 F–R 기준으로 유지
        a_dict["heterodimer_dg"] = het_fr_dg
        a_dict["heterodimer_tm"] = het_fr_tm

        # 추가: 명시적으로 세 쌍 모두 저장
        a_dict["heterodimer_fr_dg"] = het_fr_dg
        a_dict["heterodimer_fr_tm"] = het_fr_tm
        a_dict["heterodimer_fp_dg"] = het_fp_dg
        a_dict["heterodimer_fp_tm"] = het_fp_tm
        a_dict["heterodimer_rp_dg"] = het_rp_dg
        a_dict["heterodimer_rp_tm"] = het_rp_tm

        # ---------- 2) hairpin / homodimer / (F–R) heterodimer QC ----------
        hairpin_ok, homodimer_ok, hetero_fr_ok = _qc_bool_flags(a_dict, qc_thresholds)

        # ---------- 3) F–P / R–P heterodimer QC ----------
        # probe가 없는 경우에는 FP/RP QC는 "해당 없음"으로 보고, QC_PASS 계산에는 영향 없도록 True 처리.
        if f_seq and p_seq:
            hetero_fp_ok = _hetero_ok(het_fp_dg, het_fp_tm, qc_thresholds)
        else:
            hetero_fp_ok = True  # probe 없으면 이 조건은 pass 취급

        if r_seq and p_seq:
            hetero_rp_ok = _hetero_ok(het_rp_dg, het_rp_tm, qc_thresholds)
        else:
            hetero_rp_ok = True

        # ---------- 4) QC 플래그(O/X) 및 최종 QC_PASS ----------
        a_dict["QC_HAIRPIN"] = "O" if hairpin_ok else "X"
        a_dict["QC_HOMODIMER"] = "O" if homodimer_ok else "X"
        a_dict["QC_HETERODIMER_FR"] = "O" if hetero_fr_ok else "X"

        # FP/RP는 probe 유무에 따라 "O"/"X"/"-" 로 표현
        if f_seq and p_seq:
            a_dict["QC_HETERODIMER_FP"] = "O" if hetero_fp_ok else "X"
        else:
            a_dict["QC_HETERODIMER_FP"] = "-"  # probe 또는 forward 미존재

        if r_seq and p_seq:
            a_dict["QC_HETERODIMER_RP"] = "O" if hetero_rp_ok else "X"
        else:
            a_dict["QC_HETERODIMER_RP"] = "-"  # probe 또는 reverse 미존재

        qc_pass = (
            hairpin_ok
            and homodimer_ok
            and hetero_fr_ok
            and hetero_fp_ok
            and hetero_rp_ok
        )
        a_dict["QC_PASS"] = "O" if qc_pass else "X"

        total_rows.append(a_dict)

        if qc_pass:
            filtered_rows.append(a_dict)

    return total_rows, filtered_rows
