from __future__ import annotations
from typing import Any, Optional, Dict
import pysam

from pcr.seq.fetch import GenomicRegion, fetch_template_sequence
from pcr.designers.base import PrimerDesigner, ProbePrimerDesigner
from pcr.pipelines.base import run_pipeline, PipelineResult
from pcr.qc import QCThresholds

def run_qpcr(
    *,
    region: GenomicRegion,
    fasta: pysam.FastaFile,
    pcr_cfg: Any,
    qc_cfg: Any,
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

    qc_th = QCThresholds(
        hairpin_min_dg=qc_cfg.HAIRPIN_MIN_DG,
        homodimer_min_dg=qc_cfg.HOMODIMER_MIN_DG,
        heterodimer_min_dg=qc_cfg.HETERODIMER_MIN_DG,
    )
    genomic_id = f"{region.chrom}:{region.start}-{region.end}"
    return run_pipeline(genomic_id=genomic_id, designer=designer, qc_thresholds=qc_th)
