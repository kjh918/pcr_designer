# primer/core.py

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import pandas as pd
import primer3
import pysam

from app.schemas import RegionInput
from primer._get_cfg import get_fasta_handle, get_pcr_params_with_override
from config.settings import settings
from primer.designer import PrimerDesigner, ProbePrimerDesigner
from primer.qc import QCThresholds, evaluate_amplicons

# ----------------------------------------------------------------------
# 메인 엔트리: qPCR 디자인
# ----------------------------------------------------------------------
def design_qpcr_for_region(
    region: RegionInput,
    reference_name: str,
    min_amplicon_length: Optional[int] = None,
    max_amplicon_length: Optional[int] = None,
    n_probes: Optional[int] = None,
    n_primers: Optional[int] = None,
    bisulfite: Optional[bool] = None,
    # --- Primer 상세 ---
    primer_opt_length: Optional[int] = None,
    primer_min_length: Optional[int] = None,
    primer_max_length: Optional[int] = None,
    primer_opt_gc: Optional[float] = None,
    primer_min_gc: Optional[float] = None,
    primer_max_gc: Optional[float] = None,
    # --- Probe 상세 ---
    min_primer_probe_tm_diff: Optional[float] = None,
    max_primer_probe_tm_diff: Optional[float] = None,
    probe_opt_length: Optional[int] = None,
    probe_min_length: Optional[int] = None,
    probe_max_length: Optional[int] = None,
    probe_opt_tm: Optional[float] = None,
    probe_min_tm: Optional[float] = None,
    probe_max_tm: Optional[float] = None,
    probe_opt_gc: Optional[float] = None,
    probe_min_gc: Optional[float] = None,
    probe_max_gc: Optional[float] = None,
    # --- 결과 상위 n개 (필요 시)
    n_best: int = 10,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    단일 유전체 구간(region)에 대해 qPCR용 primer/probe를 설계하고,
    전체 / QC 통과된 amplicon 정보를 DataFrame으로 반환.

    반환되는 DataFrame에는:
      - Amplicon.to_dict()에서 나온 모든 필드
      - forward/reverse heterodimer_dg, heterodimer_tm
    가 포함됩니다.
    """

    # 1) high-level 파라미터: 기존 설정 + override
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

    # 2) settings 에서 primer/probe 기본 kwargs 가져오고, 웹에서 들어온 값으로 override
    pcr_cfg = settings.pcr_params

    # -------------------------
    # Primer 기본 설정
    # -------------------------
    pk = pcr_cfg.primer_kwargs

    primer_opt_length_eff = primer_opt_length if primer_opt_length is not None else pk.opt_length
    primer_min_length_eff = primer_min_length if primer_min_length is not None else pk.min_length
    primer_max_length_eff = primer_max_length if primer_max_length is not None else pk.max_length

    primer_opt_gc_eff = primer_opt_gc if primer_opt_gc is not None else pk.opt_gc
    primer_min_gc_eff = primer_min_gc if primer_min_gc is not None else pk.min_gc
    primer_max_gc_eff = primer_max_gc if primer_max_gc is not None else pk.max_gc

    primer3_global_args = dict(pk.primer3_global_args or {})

    # primer Tm 기준값 (diff 계산용, override 안 함)
    primer_base_opt_tm = primer3_global_args.get("PRIMER_OPT_TM", 60.0)

    # -------------------------
    # Probe 기본 설정 + Tm 계산
    # -------------------------
    prk = pcr_cfg.probe_kwargs

    # probe length / GC 기본값
    probe_opt_length_eff = probe_opt_length if probe_opt_length is not None else prk.opt_length
    probe_min_length_eff = probe_min_length if probe_min_length is not None else prk.min_length
    probe_max_length_eff = probe_max_length if probe_max_length is not None else prk.max_length

    probe_opt_gc_eff = probe_opt_gc if probe_opt_gc is not None else prk.opt_gc
    probe_min_gc_eff = probe_min_gc if probe_min_gc is not None else prk.min_gc
    probe_max_gc_eff = probe_max_gc if probe_max_gc is not None else prk.max_gc

    # probe Tm 기본값 (config)
    default_probe_opt_tm = prk.opt_tm
    default_probe_min_tm = prk.min_tm
    default_probe_max_tm = prk.max_tm

    # diff가 들어왔으면 primer Tm 기준으로 probe Tm 범위 계산
    if min_primer_probe_tm_diff is not None or max_primer_probe_tm_diff is not None:
        min_diff = (
            min_primer_probe_tm_diff
            if min_primer_probe_tm_diff is not None
            else max_primer_probe_tm_diff
        )
        max_diff = (
            max_primer_probe_tm_diff
            if max_primer_probe_tm_diff is not None
            else min_primer_probe_tm_diff
        )

        probe_min_tm_eff = (
            probe_min_tm
            if probe_min_tm is not None
            else primer_base_opt_tm + min_diff
        )
        probe_max_tm_eff = (
            probe_max_tm
            if probe_max_tm is not None
            else primer_base_opt_tm + max_diff
        )
        probe_opt_tm_eff = (
            probe_opt_tm
            if probe_opt_tm is not None
            else primer_base_opt_tm + (min_diff + max_diff) / 2.0
        )
    else:
        # diff 없이 절대값 또는 config 그대로
        probe_opt_tm_eff = probe_opt_tm if probe_opt_tm is not None else default_probe_opt_tm
        probe_min_tm_eff = probe_min_tm if probe_min_tm is not None else default_probe_min_tm
        probe_max_tm_eff = probe_max_tm if probe_max_tm is not None else default_probe_max_tm

    probe_primer3_global_args = dict(prk.primer3_global_args or {})

    # 3) 레퍼런스 FASTA 핸들 & 템플릿 서열 추출
    fasta = pysam.FastaFile(get_fasta_handle(reference_name))

    target_genomic_start = region.start
    target_genomic_end = region.end
    target_len = target_genomic_end - target_genomic_start + 1

    # max_amplicon_length 기준으로 타겟 앞/뒤 여유(bp) 계산
    if max_amplicon_length is not None and max_amplicon_length > target_len:
        extra = max_amplicon_length - target_len   # 타겟 밖으로 나갈 수 있는 총 길이
        upstream_extra = extra // 2
        downstream_extra = extra - upstream_extra
    else:
        # max_amplicon_length가 없거나 타겟보다 작은 경우 → 여유 없이 타겟만 사용
        upstream_extra = 0
        downstream_extra = 0

    # template 구간 (1-based, inclusive)
    template_start = max(1, target_genomic_start - upstream_extra)
    template_end = target_genomic_end + downstream_extra

    # pysam.fetch: 0-based, end-exclusive
    template_sequence: str = fasta.fetch(
        region.chrom,
        template_start - 1,
        template_end,
    ).upper()

    # 템플릿 내에서 타겟의 index (0-based, inclusive)
    target_start_index = target_genomic_start - template_start
    target_end_index = target_start_index + target_len - 1

    # reference_template_sequence 는 일반적으로 원본 서열 전체
    # (bisulfite 모드에서는 여기서 변환 전/후를 나눌 수도 있음)
    reference_template_sequence = template_sequence

    # 4) PrimerDesigner / ProbePrimerDesigner 선택해서 실행
    if n_probes is not None and n_probes > 0:
        # probe 포함 설계
        designer = ProbePrimerDesigner(
            template_sequence=template_sequence,
            target_start_index=target_start_index,
            target_end_index=target_end_index,
            reference_template_sequence=reference_template_sequence,
            min_amplicon_length=min_amplicon_length,
            max_amplicon_length=max_amplicon_length,
            n_primers=n_primers,
            n_probes=n_probes,
            # primer 조건
            opt_length=primer_opt_length_eff,
            min_length=primer_min_length_eff,
            max_length=primer_max_length_eff,
            opt_gc=primer_opt_gc_eff,
            min_gc=primer_min_gc_eff,
            max_gc=primer_max_gc_eff,
            primer3_global_args=primer3_global_args,
            # probe 조건
            probe_opt_length=probe_opt_length_eff,
            probe_min_length=probe_min_length_eff,
            probe_max_length=probe_max_length_eff,
            probe_opt_tm=probe_opt_tm_eff,
            probe_min_tm=probe_min_tm_eff,
            probe_max_tm=probe_max_tm_eff,
            probe_opt_gc=probe_opt_gc_eff,
            probe_min_gc=probe_min_gc_eff,
            probe_max_gc=probe_max_gc_eff,
            probe_primer3_global_args=probe_primer3_global_args,
        )
    else:
        # 프라이머만 설계
        designer = PrimerDesigner(
            template_sequence=template_sequence,
            target_start_index=target_start_index,
            target_end_index=target_end_index,
            reference_template_sequence=reference_template_sequence,
            min_amplicon_length=min_amplicon_length,
            max_amplicon_length=max_amplicon_length,
            n_primers=n_primers,
            opt_length=primer_opt_length_eff,
            min_length=primer_min_length_eff,
            max_length=primer_max_length_eff,
            opt_gc=primer_opt_gc_eff,
            min_gc=primer_min_gc_eff,
            max_gc=primer_max_gc_eff,
            primer3_global_args=primer3_global_args,
        )

    # PrimerDesigner / ProbePrimerDesigner 내부에서 primer3 실행
    # 메서드 이름이 design_primer() 라면 여기만 바꾸면 됨
    designer.design()

    # 5) QC Threshold 설정
    qc_th = QCThresholds(
        hairpin_max_tm=47.0,
        hairpin_min_dg=-5.0,
        homodimer_min_dg=-6.0,
        heterodimer_min_dg=-6.0,
        heterodimer_max_tm=45.0,
    )

    total_rows, filtered_rows = evaluate_amplicons(
        designer.amplicon_list,
        qc_thresholds=qc_th,
    )

    total_df = pd.DataFrame(total_rows)
    filtered_df = pd.DataFrame(filtered_rows)

    print("QC 통과 primer 개수:", len(filtered_df))
    if not filtered_df.empty:
        print(filtered_df.head())

    return total_df, filtered_df

