# pcr/config/runtime.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple, Mapping

import pysam

from pcr.config.settings import settings  # ✅ YAML 기반으로 로드된 전역 default
from pcr.config.schema.qc import QCParams
from pcr.config.schema.pcr import PCRParams  # ✅ 네 프로젝트 경로 기준 (schema/pcr.py)


# -------------------------
# settings getter (router에서 함수로 쓰고 싶을 때)
# -------------------------
def get_settings():
    return settings


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
# dict deep merge (nested override용)
# -------------------------
def deep_merge(base: Dict[str, Any], override: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
    """
    base를 복사한 새 dict에 override를 재귀 merge.
    None 값은 무시.
    """
    out: Dict[str, Any] = dict(base)
    if not override:
        return out

    for k, v in override.items():
        if v is None:
            continue
        if isinstance(v, Mapping) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


# -------------------------
# QC params merge (settings + web override)
# -------------------------
_QC_WEB_OVERRIDABLE = {
    # thermo
    "HAIRPIN_MIN_DG",
    "HOMODIMER_MIN_DG",
    "HETERODIMER_MIN_DG",
    # diff tm
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
}


def build_qc_params_from_web(overrides: Optional[Dict[str, Any]] = None) -> QCParams:
    """
    settings.qc_params + web overrides => 최종 QCParams

    - settings는 절대 수정하지 않음
    - 허용된 키만 반영
    """
    base: Dict[str, Any] = settings.qc_params.model_dump()

    if overrides:
        clean: Dict[str, Any] = {}
        for k, v in overrides.items():
            if v is None:
                continue
            if k in _QC_WEB_OVERRIDABLE:
                clean[k] = v
        base = deep_merge(base, clean)

    return QCParams.model_validate(base)


# -------------------------
# PCR params merge (settings + web override)
# -------------------------
# ✅ PCRParams 스키마가 nested일 가능성이 높아서 "nested override"를 받아 deep_merge 함
# ✅ 안전하게 하려면 allowlist를 두는 게 좋음. (일단 자주 쓰는 항목만)
_PRIMER3_WEB_OVERRIDABLE = {
    # primer
    "PRIMER_OPT_SIZE", "PRIMER_MIN_SIZE", "PRIMER_MAX_SIZE",
    "PRIMER_OPT_TM", "PRIMER_MIN_TM", "PRIMER_MAX_TM",
    "PRIMER_OPT_GC_PERCENT", "PRIMER_MIN_GC", "PRIMER_MAX_GC",
    # product size
    "PRIMER_PRODUCT_SIZE_RANGE",
    # probe(내부 올리고)
    "PRIMER_INTERNAL_OPT_SIZE", "PRIMER_INTERNAL_MIN_SIZE", "PRIMER_INTERNAL_MAX_SIZE",
    "PRIMER_INTERNAL_OPT_TM", "PRIMER_INTERNAL_MIN_TM", "PRIMER_INTERNAL_MAX_TM",
    "PRIMER_INTERNAL_MIN_GC", "PRIMER_INTERNAL_MAX_GC",
}

def _filter_primer3_args(d: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not d:
        return None
    out: Dict[str, Any] = {}
    for k, v in d.items():
        if v is None:
            continue
        if k in _PRIMER3_WEB_OVERRIDABLE:
            out[k] = v
    return out


def build_pcr_params_from_web(overrides: Optional[Dict[str, Any]] = None) -> PCRParams:
    """
    settings.pcr_params + web overrides => 최종 PCRParams (request-scope)

    overrides는 nested dict를 권장:
    {
      "primer_kwargs": {
         "min_amplicon_length": 80,
         "max_amplicon_length": 200,
         "n_primers": 5,
         "primer3_global_args": {...}
      },
      "probe_kwargs": {
         "n_probes": 1,
         "primer3_global_args": {...}
      },
      "bisulfite": {"run": False}
    }

    - settings는 절대 수정하지 않음
    - primer3_global_args는 allowlist 필터 적용(안전)
    """
    base: Dict[str, Any] = settings.pcr_params.model_dump()
    if not overrides:
        return PCRParams.model_validate(base)

    clean = dict(overrides)

    # primer3 args allowlist 적용(있는 경우만)
    try:
        pk = clean.get("primer_kwargs") or {}
        if isinstance(pk, dict) and "primer3_global_args" in pk:
            pk = dict(pk)
            pk["primer3_global_args"] = _filter_primer3_args(pk.get("primer3_global_args")) or {}
            clean["primer_kwargs"] = pk

        prk = clean.get("probe_kwargs") or {}
        if isinstance(prk, dict) and "primer3_global_args" in prk:
            prk = dict(prk)
            prk["primer3_global_args"] = _filter_primer3_args(prk.get("primer3_global_args")) or {}
            clean["probe_kwargs"] = prk
    except Exception:
        # 스키마가 다르더라도 deep_merge + model_validate에서 잡히게 둠
        pass

    merged = deep_merge(base, clean)
    return PCRParams.model_validate(merged)


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
    pcr_cfg: PCRParams,
    min_amplicon_length: Optional[int] = None,
    max_amplicon_length: Optional[int] = None,
    n_probes: Optional[int] = None,
    n_primers: Optional[int] = None,
    bisulfite: Optional[bool] = None,
) -> PCRResolvedParams:
    """
    ✅ IMPORTANT:
    - settings를 직접 보지 않고, 반드시 request-scope인 pcr_cfg를 기준으로 resolve
    """
    return PCRResolvedParams(
        min_amplicon_length=(
            min_amplicon_length
            if min_amplicon_length is not None
            else pcr_cfg.primer_kwargs.min_amplicon_length
        ),
        max_amplicon_length=(
            max_amplicon_length
            if max_amplicon_length is not None
            else pcr_cfg.primer_kwargs.max_amplicon_length
        ),
        n_primers=(n_primers if n_primers is not None else pcr_cfg.primer_kwargs.n_primers),
        n_probes=(n_probes if n_probes is not None else pcr_cfg.probe_kwargs.n_probes),
        bisulfite=(bisulfite if bisulfite is not None else pcr_cfg.bisulfite.run),
    )


def merge_dict(base: Optional[Dict[str, Any]], override: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """primer3_global_args 같은 dict merge용 (단순 merge)"""
    out: Dict[str, Any] = dict(base or {})
    if override:
        out.update({k: v for k, v in override.items() if v is not None})
    return out
