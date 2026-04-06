"""
pcr/config/loader.py
base_pcr.yaml (분석 설정)과 system.yaml (인프라 경로 설정)을 
각각 로드하여 하나의 PipelineConfig 객체로 병합합니다.
"""
import yaml
import collections.abc
import os
from typing import Dict, Any, Optional

# 스키마 경로에 맞게 import (사용자님의 실제 구조에 맞게 조절 가능)
from .schema.app import PipelineConfig

def deep_update(d: dict, u: dict) -> dict:
    """중첩 딕셔너리를 재귀적으로 병합합니다."""
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping) and isinstance(d.get(k), collections.abc.Mapping):
            d[k] = deep_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def load_pipeline_config(
    base_yaml_path: str,       # ★ 여기가 수정된 포인트입니다!
    system_yaml_path: str,     # ★ 시스템 경로 설정 파일 추가
    assay_type: str = "base", 
    user_overrides: Optional[Dict[str, Any]] = None
) -> PipelineConfig:
    
    # 1. YAML 파일 각각 읽기
    with open(base_yaml_path, 'r', encoding='utf-8') as f:
        base_yaml = yaml.safe_load(f) or {}
        
    with open(system_yaml_path, 'r', encoding='utf-8') as f:
        system_yaml = yaml.safe_load(f) or {}

    # 2. 두 설정 파일의 베이스 파라미터 병합
    merged_config = {
        "system": system_yaml.get("system", {}),
        "references": system_yaml.get("references", {}),
        "pcr_params": base_yaml.get("pcr_params", {}),
        "qc_criteria": base_yaml.get("qc_criteria", {})
    }

    # 3. Assay Overrides 병합 (base_pcr.yaml 기준)
    assay_overrides = base_yaml.get("assay_overrides", {})
    if assay_type in assay_overrides:
        specific_overrides = assay_overrides[assay_type]
        merged_config = deep_update(merged_config, specific_overrides)

    # 4. User Overrides (CLI/Web 입력) 병합
    if user_overrides:
        merged_config = deep_update(merged_config, user_overrides)

    # 5. Pydantic 타입 검증 및 객체 생성
    # (Pydantic V2 문법에 맞춰 PipelineConfig 모델링)
    return PipelineConfig(**merged_config)