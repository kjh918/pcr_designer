# primer/qc/thermo.py
from __future__ import annotations

from typing import Dict, Tuple, List, Iterable, Any

import primer3
from Bio.Seq import Seq

from pcr.qc.config import QCThresholds


# -----------------------------
# Thermo util
# -----------------------------
def compute_heterodimer(f_seq: str, r_seq: str) -> Tuple[float, float]:
    hetero = primer3.calc_heterodimer(f_seq, r_seq)
    if hetero.structure_found:
        het_dg = hetero.dg / 1000.0
        het_tm = hetero.tm
    else:
        het_dg = 0.0
        het_tm = 0.0
    return het_dg, het_tm


def compute_hairpin(seq: str) -> Tuple[float, float]:
    hp = primer3.calc_hairpin(seq)
    if hp.structure_found:
        hp_dg = hp.dg / 1000.0
        hp_tm = hp.tm
    else:
        hp_dg = 0.0
        hp_tm = 0.0
    return hp_tm, hp_dg


def compute_homodimer(seq: str) -> float:
    hd = primer3.calc_homodimer(seq)
    if hd.structure_found:
        hd_dg = hd.dg / 1000.0
    else:
        hd_dg = 0.0
    return hd_dg


def _qc_bool_flags(amp: Dict[str, Any], th: QCThresholds) -> Tuple[bool, bool, bool]:
    hairpin_ok = (
        amp.get("forward_hairpin_dg", 0.0) >= th.hairpin_min_dg
        and amp.get("reverse_hairpin_dg", 0.0) >= th.hairpin_min_dg
    )

    homodimer_ok = (
        amp.get("forward_homodimer_dg", 0.0) >= th.homodimer_min_dg
        and amp.get("reverse_homodimer_dg", 0.0) >= th.homodimer_min_dg
    )

    hetero_fr_ok = (amp.get("heterodimer_dg", 0.0) >= th.heterodimer_min_dg)
    return hairpin_ok, homodimer_ok, hetero_fr_ok


def _hetero_ok(dg: float, tm: float, th: QCThresholds) -> bool:
    return dg >= th.heterodimer_min_dg


def amplicon_passes_qc(amp: Dict[str, Any], th: QCThresholds) -> bool:
    hairpin_ok, homodimer_ok, hetero_fr_ok = _qc_bool_flags(amp, th)
    return hairpin_ok and homodimer_ok and hetero_fr_ok


def evaluate_amplicons(
    genomic_id,
    amplicons: Iterable[Any],
    qc_thresholds: QCThresholds,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    total_rows: List[Dict[str, Any]] = []
    filtered_rows: List[Dict[str, Any]] = []

    for amplicon in amplicons:
        a_dict = amplicon.to_dict()
        a_dict["ID"] = genomic_id

        f_seq = a_dict.get("forward_sequence")
        r_seq = a_dict.get("reverse_sequence")
        p_seq = a_dict.get("probe_sequence")

        # reverse RC 저장(원 코드 유지)
        if r_seq:
            a_dict["rc_reverse_sequence"] = str(Seq(r_seq).reverse_complement())
        else:
            a_dict["rc_reverse_sequence"] = ""

        # ---- hairpin/homodimer 계산(원래 SimpleAmplicon.to_dict에서 하던 것 포함) ----
        if f_seq:
            f_hp_tm, f_hp_dg = compute_hairpin(f_seq)
            f_hd_dg = compute_homodimer(f_seq)
        else:
            f_hp_tm, f_hp_dg, f_hd_dg = 0.0, 0.0, 0.0

        if r_seq:
            r_hp_tm, r_hp_dg = compute_hairpin(r_seq)
            r_hd_dg = compute_homodimer(r_seq)
        else:
            r_hp_tm, r_hp_dg, r_hd_dg = 0.0, 0.0, 0.0

        if p_seq:
            p_hd_dg = compute_homodimer(p_seq)
        else:
            p_hd_dg = 0.0

        a_dict.update(
            {
                "forward_hairpin_tm": f_hp_tm,
                "forward_hairpin_dg": f_hp_dg,
                "reverse_hairpin_tm": r_hp_tm,
                "reverse_hairpin_dg": r_hp_dg,
                "forward_homodimer_dg": f_hd_dg,
                "reverse_homodimer_dg": r_hd_dg,
                "probe_homodimer_dg": p_hd_dg,
            }
        )

        # ---------- heterodimer 계산 ----------
        if f_seq and r_seq:
            het_fr_dg, het_fr_tm = compute_heterodimer(f_seq, r_seq)
        else:
            het_fr_dg, het_fr_tm = 0.0, 0.0

        if f_seq and p_seq:
            het_fp_dg, het_fp_tm = compute_heterodimer(f_seq, p_seq)
        else:
            het_fp_dg, het_fp_tm = 0.0, 0.0

        if r_seq and p_seq:
            het_rp_dg, het_rp_tm = compute_heterodimer(r_seq, p_seq)
        else:
            het_rp_dg, het_rp_tm = 0.0, 0.0

        a_dict["heterodimer_dg"] = het_fr_dg
        a_dict["heterodimer_tm"] = het_fr_tm

        a_dict["heterodimer_fr_dg"] = het_fr_dg
        a_dict["heterodimer_fr_tm"] = het_fr_tm
        a_dict["heterodimer_fp_dg"] = het_fp_dg
        a_dict["heterodimer_fp_tm"] = het_fp_tm
        a_dict["heterodimer_rp_dg"] = het_rp_dg
        a_dict["heterodimer_rp_tm"] = het_rp_tm

        # ---------- QC flags ----------
        hairpin_ok, homodimer_ok, hetero_fr_ok = _qc_bool_flags(a_dict, qc_thresholds)

        hetero_fp_ok = _hetero_ok(het_fp_dg, het_fp_tm, qc_thresholds) if (f_seq and p_seq) else True
        hetero_rp_ok = _hetero_ok(het_rp_dg, het_rp_tm, qc_thresholds) if (r_seq and p_seq) else True

        a_dict["QC_HAIRPIN"] = "O" if hairpin_ok else "X"
        a_dict["QC_HOMODIMER"] = "O" if homodimer_ok else "X"
        a_dict["QC_HETERODIMER_FR"] = "O" if hetero_fr_ok else "X"
        a_dict["QC_HETERODIMER_FP"] = ("O" if hetero_fp_ok else "X") if (f_seq and p_seq) else "-"
        a_dict["QC_HETERODIMER_RP"] = ("O" if hetero_rp_ok else "X") if (r_seq and p_seq) else "-"

        qc_pass = hairpin_ok and homodimer_ok and hetero_fr_ok and hetero_fp_ok and hetero_rp_ok
        a_dict["QC_PASS"] = "O" if qc_pass else "X"

        total_rows.append(a_dict)
        if qc_pass:
            filtered_rows.append(a_dict)

    return total_rows, filtered_rows
