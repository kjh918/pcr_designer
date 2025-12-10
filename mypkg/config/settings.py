from __future__ import annotations

from pathlib import Path
from typing import Dict, Any

import yaml
from pydantic import BaseModel, Field


# -------------------------
# Reference 설정
# -------------------------
class ReferenceConfig(BaseModel):
    name: str
    fasta: Path
    blast: Path


# -------------------------
# PCR 관련 설정
# -------------------------
class PrimerKwargs(BaseModel):
    min_amplicon_length: int
    max_amplicon_length: int
    n_primers: int

    opt_length: int
    min_length: int
    max_length: int

    opt_gc: float
    min_gc: float
    max_gc: float

    # 비어 있을 수도 있으니 기본값 {} 허용 (mutable 주의 → default_factory)
    primer3_global_args: Dict[str, Any] = Field(default_factory=dict)


class ProbeKwargs(BaseModel):
    n_probes: int
    min_primer_probe_tm_diff: float
    max_primer_probe_tm_diff: float

    opt_length: int
    min_length: int
    max_length: int

    opt_tm: float
    min_tm: float
    max_tm: float

    opt_gc: float
    min_gc: float
    max_gc: float

    # 비어 있을 수도 있으니 기본값 {} 허용
    primer3_global_args: Dict[str, Any] = Field(default_factory=dict)


class BisulfiteConfig(BaseModel):
    run: bool
    cpg_default: str


class PCRParams(BaseModel):
    primer_kwargs: PrimerKwargs
    probe_kwargs: ProbeKwargs
    bisulfite: BisulfiteConfig


# -------------------------
# QC / BLAST 관련 설정
# -------------------------
class QCParams(BaseModel):
    # Thermo / dimer 관련 threshold
    MAX_TM_DIFF: float
    HAIRPIN_MIN_DG: float
    HOMODIMER_MIN_DG: float
    HETERODIMER_MIN_DG: float

    # BLAST 및 amplicon 관련 설정
    BLAST_ROOT: Path
    BLAST_BIN_DIR: Path
    BLASTN: Path
    BLASTDBCMD: Path

    BLAST_IDENTITY_THRESHOLD: float
    BLAST_MAX_ALIGNMENTS: int
    BLAST_LENGTH_THRESHOLD: int
    MIN_AMP_BP: int
    MAX_AMP_BP: int


# -------------------------
# 전체 Settings
# -------------------------
class Settings(BaseModel):
    references: Dict[str, ReferenceConfig]
    pcr_params: PCRParams
    qc_params: QCParams


# -------------------------
# YAML 로드 함수
# -------------------------
def load_settings() -> Settings:
    """
    프로젝트 루트(BASE_DIR)/config/parameter.yaml 을 로드해서
    Settings 객체로 반환.
    """
    BASE_DIR = Path(__file__).resolve().parents[1]
    CONFIG_FILE = BASE_DIR / "config" / "parameter.yaml"

    if not CONFIG_FILE.exists():
        raise FileNotFoundError(f"Config file not found: {CONFIG_FILE}")

    with open(CONFIG_FILE, "r") as f:
        raw = yaml.safe_load(f) or {}

    # ---- references ----
    refs: Dict[str, ReferenceConfig] = {}
    for name, cfg in (raw.get("references") or {}).items():
        refs[name] = ReferenceConfig(
            name=name,
            fasta=Path(cfg["fasta"]),
            blast=Path(cfg["blast"]),
        )

    # ---- pcr_params ----
    pcr_raw = raw.get("pcr_params") or {}
    pcr_params = PCRParams(
        primer_kwargs=PrimerKwargs(**(pcr_raw.get("primer_kwargs") or {})),
        probe_kwargs=ProbeKwargs(**(pcr_raw.get("probe_kwargs") or {})),
        bisulfite=BisulfiteConfig(**(pcr_raw.get("bisulfite") or {})),
    )

    # ---- qc_params ----
    # YAML 키 이름이 'qc_params' 또는 예전 이름 'qc_thresholds' 일 수 있다고 가정
    qc_raw = (
        raw.get("qc_params")
        or {}
    )

    # Path 필드는 여기서 한 번 Path로 변환해 주는 게 안전
    if "BLAST_ROOT" in qc_raw:
        qc_raw["BLAST_ROOT"] = Path(qc_raw["BLAST_ROOT"])
    if "BLAST_BIN_DIR" in qc_raw:
        qc_raw["BLAST_BIN_DIR"] = Path(qc_raw["BLAST_BIN_DIR"])
    if "BLASTN" in qc_raw:
        qc_raw["BLASTN"] = Path(qc_raw["BLASTN"])
    if "BLASTDBCMD" in qc_raw:
        qc_raw["BLASTDBCMD"] = Path(qc_raw["BLASTDBCMD"])

    qc_params = QCParams(**qc_raw)

    return Settings(
        references=refs,
        pcr_params=pcr_params,
        qc_params=qc_params,
    )


settings = load_settings()
