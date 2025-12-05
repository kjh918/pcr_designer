# primer/qc.py

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, List, Iterable, Any

import primer3


def compute_heterodimer(f_seq: str, r_seq: str) -> Tuple[float, float]:
    """
    Forward / Reverse 프라이머의 heterodimer ΔG / Tm 계산.
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


def amplicon_passes_qc(amp: Dict[str, float], th: QCThresholds) -> bool:
    """
    amplicon.to_dict() 결과 + heterodimer 정보가 들어있는 dict를 받아
    QC 통과 여부를 판단.
    """
    hairpin_ok = (
        amp.get("forward_hairpin_tm", 0.0) <= th.hairpin_max_tm
        and amp.get("reverse_hairpin_tm", 0.0) <= th.hairpin_max_tm
        and amp.get("forward_hairpin_dg", 0.0) >= th.hairpin_min_dg
        and amp.get("reverse_hairpin_dg", 0.0) >= th.hairpin_min_dg
    )

    homodimer_ok = (
        amp.get("forward_homodimer_dg", 0.0) >= th.homodimer_min_dg
        and amp.get("reverse_homodimer_dg", 0.0) >= th.homodimer_min_dg
    )

    heterodimer_ok = (
        amp.get("heterodimer_dg", 0.0) >= th.heterodimer_min_dg
        and amp.get("heterodimer_tm", 0.0) <= th.heterodimer_max_tm
    )

    return hairpin_ok and homodimer_ok and heterodimer_ok


def evaluate_amplicons(
    amplicons: Iterable[Any],
    qc_thresholds: QCThresholds,
) -> Tuple[List[Dict[str, float]], List[Dict[str, float]]]:
    """
    Amplicon 객체 리스트에 대해:
      - Amplicon.to_dict() 호출
      - heterodimer_dg / heterodimer_tm 계산 추가
      - QC 통과 여부 필터링

    Parameters
    ----------
    amplicons : Iterable[Amplicon-like]
        각 원소가 .to_dict() 메서드를 가지고 있다고 가정.
    qc_thresholds : QCThresholds
        QC 기준 값.

    Returns
    -------
    total_rows, filtered_rows : (List[Dict], List[Dict])
        - total_rows: 모든 amplicon의 dict
        - filtered_rows: QC 통과한 amplicon의 dict
    """
    total_rows: List[Dict[str, float]] = []
    filtered_rows: List[Dict[str, float]] = []

    for amplicon in amplicons:
        a_dict = amplicon.to_dict()

        f_seq = a_dict.get("forward_sequence")
        r_seq = a_dict.get("reverse_sequence")

        if f_seq and r_seq:
            het_dg, het_tm = compute_heterodimer(f_seq, r_seq)
        else:
            het_dg, het_tm = 0.0, 0.0

        a_dict["heterodimer_dg"] = het_dg
        a_dict["heterodimer_tm"] = het_tm

        total_rows.append(a_dict)

        if amplicon_passes_qc(a_dict, qc_thresholds):
            filtered_rows.append(a_dict)

    return total_rows, filtered_rows
