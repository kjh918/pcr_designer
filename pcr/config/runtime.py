# pcr/config/runtime.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

import pysam

from pcr.config.settings import settings  # ✅ 이미 로드된 settings를 재사용
from pcr.config.schema.qc import QCParams


# -------------------------
# FASTA handle cache (pysam only)
# -------------------------
_FASTA_CACHE: Dict[str, pysam.FastaFile] = {}


def get_reference_paths(ref_name: str) -> Tuple[str, str]:
    """(fasta_path, blast_db_path)"""
    if ref_name not in settings.references:
        raise ValueError(f"Unknown reference name: {ref_name}")
    ref = settings.references[ref_name]
    return str(ref.fasta), str(ref.blast)


def get_fasta_handle(ref_name: str) -> pysam.FastaFile:
    """ref_name -> pysam.FastaFile (cached)"""
    if ref_name not in _FASTA_CACHE:
        fasta_path, _ = get_reference_paths(ref_name)
        _FASTA_CACHE[ref_name] = pysam.FastaFile(fasta_path)
    return _FASTA_CACHE[ref_name]


def clear_fasta_cache(ref_name: Optional[str] = None) -> None:
    """테스트/운영에서 캐시 초기화 필요할 때"""
    if ref_name is None:
        _FASTA_CACHE.clear()
    else:
        _FASTA_CACHE.pop(ref_name, None)


# -------------------------
# QC params merge (settings + web override)
# -------------------------
_QC_WEB_OVERRIDABLE = {
    # thermo
    "HAIRPIN_MIN_DG",
    "HOMODIMER_MIN_DG",
    "HETERODIMER_MIN_DG",
    # diff tm (템플릿에서 조정 가능하게)
    "PRIMER_MAX_DIFF_TM",
    "PRIMER_MIN_DIFF_TM",
    "PROBE_MAX_DIFF_TM",
    "PROBE_MIN_DIFF_TM",
    # blast filter
    "BLAST_IDENTITY_THRESHOLD",
    "BLAST_LENGTH_THRESHOLD",
    "BLAST_MAX_ALIGNMENTS",
    "MIN_AMP_BP",
    "MAX_AMP_BP",
    # 보통 실행 경로는 웹에서 안 바꾸는 게 안정적이라 제외
    # 필요하면 아래도 추가 가능:
    # "BLAST_ROOT", "BLAST_BIN_DIR", "BLASTN", "BLASTDBCMD",
}


def build_qc_params_from_web(overrides: Optional[Dict[str, Any]] = None) -> QCParams:
    """
    settings.qc_params + web overrides => 최종 QCParams

    주의:
    - pcr/qc 내부에서는 settings import 금지
    - router/pipeline에서 이 함수로 qc_params를 만들어 주입
    """
    base: Dict[str, Any] = settings.qc_params.model_dump()

    if overrides:
        for k, v in overrides.items():
            if v is None:
                continue
            if k in _QC_WEB_OVERRIDABLE:
                base[k] = v

    return QCParams.model_validate(base)


# -------------------------
# PCR param resolve helpers
# -------------------------
@dataclass(frozen=True)
class PCRResolvedParams:
    min_amplicon_length: int
    max_amplicon_length: int
    n_probes: int
    n_primers: int
    bisulfite: bool


def resolve_pcr_params(
    *,
    min_amplicon_length: Optional[int] = None,
    max_amplicon_length: Optional[int] = None,
    n_probes: Optional[int] = None,
    n_primers: Optional[int] = None,
    bisulfite: Optional[bool] = None,
) -> PCRResolvedParams:
    """
    override(None 허용) + settings 기본값으로 최종 파라미터 확정
    """
    p = settings.pcr_params
    return PCRResolvedParams(
        min_amplicon_length=(
            min_amplicon_length
            if min_amplicon_length is not None
            else p.primer_kwargs.min_amplicon_length
        ),
        max_amplicon_length=(
            max_amplicon_length
            if max_amplicon_length is not None
            else p.primer_kwargs.max_amplicon_length
        ),
        n_primers=(n_primers if n_primers is not None else p.primer_kwargs.n_primers),
        n_probes=(n_probes if n_probes is not None else p.probe_kwargs.n_probes),
        bisulfite=(bisulfite if bisulfite is not None else p.bisulfite.run),
    )


def merge_dict(base: Optional[Dict[str, Any]], override: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """primer3_global_args 같은 dict merge용"""
    out: Dict[str, Any] = dict(base or {})
    if override:
        out.update(override)
    return out
