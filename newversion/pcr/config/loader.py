import yaml
import os
from typing import Dict, Any, Optional
from .schema.root import AppConfig

class ConfigLoader:
    def __init__(self, config_dir: Optional[str] = None):
        # 현재 loader.py 파일의 위치를 기준으로 경로 설정
        if config_dir is None:
            self.config_dir = os.path.dirname(os.path.abspath(__file__))
        else:
            self.config_dir = config_dir

    def _load_yaml(self, path: str) -> Dict[str, Any]:
        """YAML 파일을 읽어 딕셔너리로 반환"""
        if not os.path.exists(path):
            raise FileNotFoundError(f"설정 파일을 찾을 수 없습니다: {path}")
        with open(path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f) or {}

    def _deep_merge(self, base: dict, override: dict) -> dict:
        """재귀적으로 딕셔너리를 병합 (Deep Merge)"""
        for key, value in override.items():
            if isinstance(value, dict) and key in base and isinstance(base[key], dict):
                self._deep_merge(base[key], value)
            else:
                base[key] = value
        return base

    def load(self, preset_name: str = "default") -> AppConfig:
        """
        1. system.yaml (경로, 레퍼런스) 로드
        2. presets/default.yaml 로드
        3. 선택된 preset (예: high_gc.yaml)으로 덮어쓰기
        4. AppConfig 객체로 변환
        """
        # 1. 시스템 공통 설정
        system_path = os.path.join(self.config_dir, "system.yaml")
        final_dict = self._load_yaml(system_path)

        # 2. 기본 프리셋 로드
        default_preset_path = os.path.join(self.config_dir, "presets", "default.yaml")
        default_dict = self._load_yaml(default_preset_path)

        # 3. 프리셋 병합 로직
        if preset_name != "default":
            preset_path = os.path.join(self.config_dir, "presets", f"{preset_name}.yaml")
            override_dict = self._load_yaml(preset_path)
            # default 위에 특정 preset 내용을 덮어씀
            default_dict = self._deep_merge(default_dict, override_dict)

        # 4. 시스템 설정과 프리셋 설정 합치기
        # final_dict(system.yaml)에 default_dict 내용을 합산
        final_dict = self._deep_merge(final_dict, default_dict)

        # 5. Pydantic 검증 및 객체화
        return AppConfig(**final_dict)

def load_config(preset_name: str = "default") -> AppConfig:
    """외부에서 간편하게 호출하기 위한 헬퍼 함수"""
    return ConfigLoader().load(preset_name)