# primer/core.py

from typing import Optional, Tuple
from app.schemas import RegionInput, PrimerPair, Primer
from primer.qpcr_designer import qPCRdesigner
from primer._get_cfg import get_fasta_handle, get_pcr_params_with_override
from config.settings import settings

import pysam


def _primer_obj_to_schema(primer_obj) -> Optional[Primer]:
    """
    pcr_components.Primer / Amplicon 내 primer 객체를
    app.schemas.Primer 로 변환.
    """
    if primer_obj is None:
        return None

    if hasattr(primer_obj, "sequence"):
        seq = primer_obj.sequence
    elif hasattr(primer_obj, "seq"):
        seq = primer_obj.seq
    else:
        raise ValueError("Primer object must have 'sequence' or 'seq' attribute")

    tm = getattr(primer_obj, "tm", None)
    gc = getattr(primer_obj, "gc", None)
    start = getattr(primer_obj, "start", None)
    end = getattr(primer_obj, "end", None)
    strand = getattr(primer_obj, "strand", None)

    return Primer(
        seq=seq,
        tm=tm,
        gc=gc,
        start=start,
        end=end,
        strand=strand,
    )

def design_qpcr_for_region(
        region: RegionInput,
        reference_name: str,
        min_amplicon_length: int | None = None,
        max_amplicon_length: int | None = None,
        n_probes: int | None = None,
        n_primers: int | None = None,
        bisulfite: bool | None = None,

        # --- Primer 상세 ---
        primer_opt_length: int | None = None,
        primer_min_length: int | None = None,
        primer_max_length: int | None = None,
        primer_opt_gc: float | None = None,
        primer_min_gc: float | None = None,
        primer_max_gc: float | None = None,

        # --- Probe 상세 ---
        probe_opt_length: int | None = None,
        probe_min_length: int | None = None,
        probe_max_length: int | None = None,
        probe_opt_tm: float | None = None,
        probe_min_tm: float | None = None,
        probe_max_tm: float | None = None,
        probe_opt_gc: float | None = None,
        probe_min_gc: float | None = None,
        probe_max_gc: float | None = None,

        n_best: int = 10,
    ):
    # 1) high-level (min/max amplicon, n_* , bisulfite) 는 기존처럼 get_pcr_params_with_override 사용
    (
        min_amplicon_length,
        max_amplicon_length,
        n_probes,
        n_primers,
        bisulfite,
    ) = get_pcr_params_with_override(
        min_amplicon_length=min_amplicon_length,
        max_amplicon_length=max_amplicon_length,
        n_probes=n_probes,
        n_primers=n_primers,
        bisulfite=bisulfite,
    )

    # 2) settings.pcr_params.primer_kwargs / probe_kwargs 를 가져와서
    #    웹에서 넘어온 값이 None이 아니면 그 값으로 override
    pcr_cfg = settings.pcr_params 
    
    pk = pcr_cfg.primer_kwargs
    primer_kwargs = {
        "opt_length": primer_opt_length if primer_opt_length is not None else pk.opt_length,
        "min_length": primer_min_length if primer_min_length is not None else pk.min_length,
        "max_length": primer_max_length if primer_max_length is not None else pk.max_length,
        "opt_gc": primer_opt_gc if primer_opt_gc is not None else pk.opt_gc,
        "min_gc": primer_min_gc if primer_min_gc is not None else pk.min_gc,
        "max_gc": primer_max_gc if primer_max_gc is not None else pk.max_gc,
        "primer3_global_args": pk.primer3_global_args,
    }

    prk = pcr_cfg.probe_kwargs
    probe_kwargs = {
        "n_probes": n_probes,
        "opt_length": probe_opt_length if probe_opt_length is not None else prk.opt_length,
        "min_length": probe_min_length if probe_min_length is not None else prk.min_length,
        "max_length": probe_max_length if probe_max_length is not None else prk.max_length,
        "opt_tm": probe_opt_tm if probe_opt_tm is not None else prk.opt_tm,
        "min_tm": probe_min_tm if probe_min_tm is not None else prk.min_tm,
        "max_tm": probe_max_tm if probe_max_tm is not None else prk.max_tm,
        "opt_gc": probe_opt_gc if probe_opt_gc is not None else prk.opt_gc,
        "min_gc": probe_min_gc if probe_min_gc is not None else prk.min_gc,
        "max_gc": probe_max_gc if probe_max_gc is not None else prk.max_gc,
        "primer3_global_args": prk.primer3_global_args,
    }

    fasta = pysam.FastaFile(get_fasta_handle(reference_name))
    
    qp = qPCRdesigner(
        f_reference_fasta=fasta,
        chrom=region.chrom,
        start=region.start,
        end=region.end,
        min_amplicon_length=min_amplicon_length,
        max_amplicon_length=max_amplicon_length,
        n_probes=n_probes,
        n_primers=n_primers,
        bisulfite=bisulfite,
        primer_kwargs=primer_kwargs,
        probe_kwargs=probe_kwargs,
    )
