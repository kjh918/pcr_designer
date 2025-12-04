# primer/core.py

from typing import Optional, Tuple
from app.schemas import RegionInput, PrimerPair, Primer
from primer.qpcr_designer import qPCRdesigner
from primer._get_cfg import get_fasta_handle, get_pcr_params_with_override
from config.settings import settings

import pysam
import pandas as pd 
import primer3

def compute_heterodimer(f_seq: str, r_seq: str):
    """F/R heterodimer ΔG / Tm 계산."""
    hetero = primer3.calc_heterodimer(f_seq, r_seq)
    if hetero.structure_found:
        het_dg = hetero.dg / 1000.0
        het_tm = hetero.tm
    else:
        het_dg = 0.0
        het_tm = 0.0
    return het_dg, het_tm

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

def amplicon_passes_qc(a: dict, th: dict) -> bool:
    """
    a: amplicon.to_dict() 결과 (dict)
    th: qc_thresholds dict
    """

    # 1) Primer Tm (각각 범위 + F/R ΔTm)
    # 3) Hairpin (Tm, dG 기준)
    hairpin_ok = (
        a["forward_hairpin_tm"] <= th["hairpin_max_tm"]
        and a["reverse_hairpin_tm"] <= th["hairpin_max_tm"]
        and a["forward_hairpin_dg"] >= th["hairpin_min_dg"]
        and a["reverse_hairpin_dg"] >= th["hairpin_min_dg"]
    )

    # 4) Homodimer (각각 dG 기준)
    homodimer_ok = (
        a["forward_homodimer_dg"] >= th["homodimer_min_dg"]
        and a["reverse_homodimer_dg"] >= th["homodimer_min_dg"]
    )

    # 5) Heterodimer (dG, Tm 기준)
    heterodimer_ok = (
        a["heterodimer_dg"] >= th["heterodimer_min_dg"]
        and a["heterodimer_tm"] <= th["heterodimer_max_tm"]
    )

    return hairpin_ok and homodimer_ok and heterodimer_ok


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
    qc_thresholds = {
        "hairpin_max_tm": 47.0,
        "hairpin_min_dg": -5.0,
        "homodimer_min_dg": -6.0,
        "heterodimer_min_dg": -6.0,
        "heterodimer_max_tm": 45.0,
    }

    filtered_amplicons = []
    filtered_rows = []
    total_rows = []

    for amplicon in qp.amplicon_list:
        a_dict = amplicon.to_dict()
        f_primer, r_primer = a_dict['forward_sequence'], a_dict['reverse_sequence']

        het_dg, het_tm = compute_heterodimer(f_primer, r_primer)
        a_dict['heterodimer_dg'] = het_dg 
        a_dict['heterodimer_tm'] = het_tm 
        amplicon_df = pd.DataFrame(a_dict.items()).T
        total_rows.append(a_dict)
        # print(amplicon_df)  # 디버깅용

        if amplicon_passes_qc(a_dict, qc_thresholds):
            filtered_amplicons.append(amplicon)
            filtered_rows.append(a_dict)

    filtered_df = pd.DataFrame(filtered_rows)
    total_df    = pd.DataFrame(total_rows)
    print("QC 통과 primer 개수:", len(filtered_df))
    print(filtered_df.head())
    return total_df, filtered_df