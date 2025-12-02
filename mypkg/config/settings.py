from pathlib import Path
from typing import Dict
import yaml
from pydantic import BaseModel
import os


# -------------------------
# Reference 설정
# -------------------------
class ReferenceConfig(BaseModel):
    name: str
    fasta: Path

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

    # 비어 있을 수도 있으니 기본값 {} 허용
    primer3_global_args: Dict[str, str] = {}


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

    # 주석 처리될 수 있으니 기본값 {} 허용
    primer3_global_args: Dict[str, str] = {}

class BisulfiteConfig(BaseModel):
    run: bool
    cpg_default: str

class PCRParams(BaseModel):
    primer_kwargs: PrimerKwargs
    probe_kwargs: ProbeKwargs
    bisulfite: BisulfiteConfig

# -------------------------
# 전체 Settings
# -------------------------
class Settings(BaseModel):
    references: Dict[str, ReferenceConfig]
    pcr_params: PCRParams

# -------------------------
# YAML 로드 함수
# -------------------------
def load_settings() -> Settings:
    """
    프로젝트 루트(BASE_DIR)/config/parameter.yaml 을 로드해서
    Settings 객체로 반환
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
        # fasta가 절대경로이므로 그대로 Path 처리
        refs[name] = ReferenceConfig(
            name=name,
            fasta=Path(cfg["fasta"]),
        )
        
    # ---- pcr_params ----
    pcr_raw = raw.get("pcr_params") or {}

    pcr_params = PCRParams(
        primer_kwargs=PrimerKwargs(**(pcr_raw.get("primer_kwargs") or {})),
        probe_kwargs=ProbeKwargs(**(pcr_raw.get("probe_kwargs") or {})),
        bisulfite=BisulfiteConfig(**(pcr_raw.get("bisulfite") or {})),
    )

    return Settings(
        references=refs,
        pcr_params=pcr_params,
    )


settings = load_settings()
