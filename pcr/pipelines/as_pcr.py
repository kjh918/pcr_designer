# pcr/pipelines/as_pcr.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, List

import pandas as pd

from pcr.seq.fetch import GenomicRegion, fetch_template_sequence, build_ref_alt_templates
from pcr.config.schema.qc import QCParams
from pcr.designers.as_pcr import AsPcrDesigner
from pcr.pipelines.base import run_pipeline, PipelineResult

def run_as_pcr_pipeline(
    *,
    genomic_id: str,
    region: GenomicRegion,
    fasta: pysam.FastaFile,
    ref_allele: str,
    alt_allele: str,
    pcr_cfg: Any,  # (지금은 안 쓰지만 유지 가능)
    qc_params: QCParams,  # ✅ qc_cfg → qc_params 로 명확히
    min_amplicon_length: Optional[int],
    max_amplicon_length: Optional[int],
    n_primers: Optional[int],
    n_reverse: Optional[int],
    forward_max_len: Optional[int],
    forward_min_len: Optional[int],
    primer3_global_args: Optional[Dict],
) -> AsPipelineResult:
    """
    AS-PCR 후보 생성 + Thermo QC + BLAST QC까지 수행 후
    total_df / filtered_df 반환.

    filtered_df 기준:
      - WT/ALT 모두 THERMO_PASS == O
      - WT/ALT 모두 BLAST_PASS == O
    """

    ref_alt_templates_dict = build_ref_alt_templates(
        region=region,
        fasta=fasta,
        ref_allele=ref_allele,
        alt_allele=alt_allele,
        max_amplicon_length=max_amplicon_length
    )

    designer = AsPcrDesigner(
        template_sequence=ref_alt_templates_dict["alt_template_sequence"],
        reference_template_sequence=ref_alt_templates_dict["ref_template_sequence"],

        target_start_index=ref_alt_templates_dict["target_start_index"],
        target_end_index=ref_alt_templates_dict["target_start_index"],  # SNP 1bp이면 동일
        target_index=ref_alt_templates_dict["target_start_index"],      # 앵커 기준도 동일

        ref_allele=ref_allele.upper(),
        alt_allele=alt_allele.upper(),

        min_amplicon_length=min_amplicon_length,
        max_amplicon_length=max_amplicon_length,

        n_primers=int(n_reverse) if n_reverse is not None else 50,
        primer3_global_args=None,
    )
    # ✅ QCThresholds 생성 제거 → QCParams 그대로 주입
    return run_pipeline(
        genomic_id=genomic_id,
        designer=designer,
        qc_params=qc_params,
    )
