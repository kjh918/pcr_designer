from __future__ import annotations

from typing import Any, Optional, Dict

import pysam

from pcr.seq.fetch import GenomicRegion, fetch_template_sequence
from pcr.designers.base import PrimerDesigner, ProbePrimerDesigner
from pcr.pipelines.base import run_pipeline, PipelineResult
from pcr.config.schema.qc import QCParams


def run_qpcr(
    *,
    region: GenomicRegion,
    fasta: pysam.FastaFile,
    pcr_cfg: Any,  # (지금은 안 쓰지만 유지 가능)
    qc_params: QCParams,  # ✅ qc_cfg → qc_params 로 명확히
    min_amplicon_length: Optional[int],
    max_amplicon_length: Optional[int],
    n_probes: Optional[int],
    n_primers: Optional[int],
    primer3_global_args: Dict,
    probe_primer3_global_args: Dict,
    # + primer/probe override들…
) -> PipelineResult:
    template_sequence, _, _, target_start_index, target_end_index = fetch_template_sequence(
        fasta, region, max_amplicon_length=max_amplicon_length
    )
    reference_template_sequence = template_sequence

    if n_probes and n_probes > 0:
        designer = ProbePrimerDesigner(
            template_sequence=template_sequence,
            reference_template_sequence=reference_template_sequence,
            target_start_index=target_start_index,
            target_end_index=target_end_index,
            min_amplicon_length=min_amplicon_length,
            max_amplicon_length=max_amplicon_length,
            n_primers=n_primers,
            n_probes=n_probes,
            primer3_global_args=primer3_global_args,
            probe_primer3_global_args=probe_primer3_global_args,
        )
    else:
        designer = PrimerDesigner(
            template_sequence=template_sequence,
            reference_template_sequence=reference_template_sequence,
            target_start_index=target_start_index,
            target_end_index=target_end_index,
            min_amplicon_length=min_amplicon_length,
            max_amplicon_length=max_amplicon_length,
            n_primers=n_primers,
            primer3_global_args=primer3_global_args,
        )

    genomic_id = f"{region.chrom}:{region.start}-{region.end}"

    # ✅ QCThresholds 생성 제거 → QCParams 그대로 주입
    return run_pipeline(
        genomic_id=genomic_id,
        designer=designer,
        qc_params=qc_params,
    )
