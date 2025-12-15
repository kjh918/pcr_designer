# config/runtime.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

import pysam

from config.settings import settings  # ✅ 이미 로드된 settings를 재사용


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
        min_amplicon_length=min_amplicon_length if min_amplicon_length is not None else p.primer_kwargs.min_amplicon_length,
        max_amplicon_length=max_amplicon_length if max_amplicon_length is not None else p.primer_kwargs.max_amplicon_length,
        n_primers=n_primers if n_primers is not None else p.primer_kwargs.n_primers,
        n_probes=n_probes if n_probes is not None else p.probe_kwargs.n_probes,
        bisulfite=bisulfite if bisulfite is not None else p.bisulfite.run,
    )


def merge_dict(base: Optional[Dict[str, Any]], override: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """primer3_global_args 같은 dict merge용"""
    out: Dict[str, Any] = dict(base or {})
    if override:
        out.update(override)
    return out
