# pcr/pipelines/as_pcr.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, List

import pandas as pd

from pcr.config.schema.qc import QCParams
from pcr.designers.as_pcr import AsPcrDesigner
from pcr.qc.thermo import evaluate_amplicons
from pcr.qc.blast import blast_qc_for_primer_pair
from pcr.qc.factories import make_amplicon_for_qc


@dataclass(frozen=True)
class AsPipelineResult:
    genomic_id: str
    total_df: pd.DataFrame
    filtered_df: pd.DataFrame


def _prefix_dict(d: Dict[str, Any], prefix: str) -> Dict[str, Any]:
    return {f"{prefix}{k}": v for k, v in d.items()}


def _run_reaction_qc(
    *,
    reaction_id: str,
    forward_seq: str,
    reverse_seq: str,
    template_seq: Optional[str],
    blast_db: str,
    qc_params: QCParams,
) -> Dict[str, Any]:
    """
    단일 reaction(F/R) 기준:
    - Thermo QC (primer3)
    - BLAST QC
    결과 dict 반환 (Thermo row + BLAST key들)
    """
    # ---- Thermo QC ----
    amp = make_amplicon_for_qc(
        forward_seq=forward_seq,
        reverse_seq=reverse_seq,
        probe_seq=None,
        template_seq=template_seq,
    )
    total_rows, _filtered_rows = evaluate_amplicons(
        reaction_id,
        [amp],
        qc_params=qc_params,
    )
    thermo_row = total_rows[0] if total_rows else {}
    thermo_pass = (thermo_row.get("QC_PASS") == "O")

    # ---- BLAST QC ----
    blast = blast_qc_for_primer_pair(
        f_name="FORWARD",
        f_seq=forward_seq,
        r_name="REVERSE",
        r_seq=reverse_seq,
        db=blast_db,
        qc_params=qc_params,
        probe_name=None,
        probe_seq=None,
    )
    blast_pass = (
        (blast.get("qc_blast_hit") == "O")
        and (blast.get("qc_blast_amplicon") == "O")
    )

    # ---- merge ----
    return {
        **thermo_row,
        "THERMO_PASS": "O" if thermo_pass else "X",
        "BLAST_F_HITS": blast.get("f_hits"),
        "BLAST_R_HITS": blast.get("r_hits"),
        "BLAST_NEARBY_AMP_COUNT": blast.get("nearby_count"),
        "BLAST_MIN_AMP_SIZE": blast.get("min_amplicon_size"),
        "BLAST_AMP_DETAIL": ";".join(blast.get("amplicon_details") or []),
        "QC_BLAST_HIT": blast.get("qc_blast_hit"),
        "QC_BLAST_AMP": blast.get("qc_blast_amplicon"),
        "BLAST_PASS": "O" if blast_pass else "X",
    }


def run_as_pcr_pipeline(
    *,
    genomic_id: str,
    designer: AsPcrDesigner,
    qc_params: QCParams,
    blast_db: str,                 # ✅ reference에서 resolve해서 넘겨줘야 함
    n_best: Optional[int] = None,  # 후보 수 제한(선택)
) -> AsPipelineResult:
    """
    AS-PCR 후보 생성 + Thermo QC + BLAST QC까지 수행 후
    total_df / filtered_df 반환.

    filtered_df 기준:
      - WT/ALT 모두 THERMO_PASS == O
      - WT/ALT 모두 BLAST_PASS == O
    """
    assays = designer.design()

    total_rows: List[Dict[str, Any]] = []
    filtered_rows: List[Dict[str, Any]] = []

    for idx, assay in enumerate(assays):
        assay_id = f"{genomic_id}|AS{idx+1}"

        # WT reaction
        wt_qc = _run_reaction_qc(
            reaction_id=f"{assay_id}|WT",
            forward_seq=assay.wt_forward.sequence,
            reverse_seq=assay.reverse.sequence,
            template_seq=assay.amplicon_sequence,
            blast_db=blast_db,
            qc_params=qc_params,
        )

        # ALT reaction
        alt_qc = _run_reaction_qc(
            reaction_id=f"{assay_id}|ALT",
            forward_seq=assay.alt_forward.sequence,
            reverse_seq=assay.reverse.sequence,
            template_seq=assay.amplicon_sequence,
            blast_db=blast_db,
            qc_params=qc_params,
        )

        wt_thermo_ok = (wt_qc.get("THERMO_PASS") == "O")
        alt_thermo_ok = (alt_qc.get("THERMO_PASS") == "O")
        wt_blast_ok = (wt_qc.get("BLAST_PASS") == "O")
        alt_blast_ok = (alt_qc.get("BLAST_PASS") == "O")

        assay_thermo_pass = wt_thermo_ok and alt_thermo_ok
        assay_blast_pass = wt_blast_ok and alt_blast_ok
        assay_pass = assay_thermo_pass and assay_blast_pass

        base = {
            "ID": genomic_id,
            "ASSAY_ID": assay_id,

            # assay core
            "wt_forward_sequence": assay.wt_forward.sequence,
            "alt_forward_sequence": assay.alt_forward.sequence,
            "reverse_sequence": assay.reverse.sequence,

            "amplicon_start_index": assay.amplicon_start_index,
            "amplicon_end_index": assay.amplicon_end_index,
            "amplicon_length": len(assay.amplicon_sequence),
            "amplicon_sequence": assay.amplicon_sequence,

            # assay pass flags
            "ASSAY_QC_PASS_THERMO": "O" if assay_thermo_pass else "X",
            "ASSAY_QC_PASS_BLAST": "O" if assay_blast_pass else "X",
            "ASSAY_QC_PASS": "O" if assay_pass else "X",
        }

        merged = {
            **base,
            **_prefix_dict(wt_qc, "WT_"),
            **_prefix_dict(alt_qc, "ALT_"),
        }

        total_rows.append(merged)
        if assay_pass:
            filtered_rows.append(merged)

        if n_best is not None and len(total_rows) >= n_best:
            break

    total_df = pd.DataFrame(total_rows)
    filtered_df = pd.DataFrame(filtered_rows)

    if len(total_df) > 0:
        total_df.index = [genomic_id] * len(total_df)
    if len(filtered_df) > 0:
        filtered_df.index = [genomic_id] * len(filtered_df)

    designer.reset()
    return AsPipelineResult(genomic_id=genomic_id, total_df=total_df, filtered_df=filtered_df)
