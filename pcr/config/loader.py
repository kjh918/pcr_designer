from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Mapping, Optional
import yaml

from config.schema.root import Settings

def _normalize_references(raw: Mapping[str, Any]) -> Dict[str, Any]:
    return dict(raw.get("references") or {})

def _normalize_pcr_params(raw: Mapping[str, Any]) -> Dict[str, Any]:
    pcr_raw = raw.get("pcr_params") or {}
    return {
        "primer_kwargs": pcr_raw.get("primer_kwargs") or {},
        "probe_kwargs": pcr_raw.get("probe_kwargs") or {},
        "bisulfite": pcr_raw.get("bisulfite") or {},
    }

def _normalize_qc_params(raw: Mapping[str, Any]) -> Dict[str, Any]:
    return dict(raw.get("qc_params") or {})

def load_settings(config_file: Optional[Path] = None) -> Settings:
    base_dir = Path(__file__).resolve().parents[1]
    config_path = config_file or (base_dir / "config" / "parameter.yaml")

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        raw = yaml.safe_load(f) or {}

    normalized = {
        "references": _normalize_references(raw),
        "pcr_params": _normalize_pcr_params(raw),
        "qc_params": _normalize_qc_params(raw),
    }
    return Settings.model_validate(normalized)
