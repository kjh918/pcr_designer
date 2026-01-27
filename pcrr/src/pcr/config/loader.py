import yaml
from pathlib import Path
from typing import Optional
from .schema.root import AppConfig

# 기본 parameter.yaml 위치 설정 (config 패키지 상위 디렉토리 기준)
DEFAULT_CONFIG_PATH = Path(__file__).parent.parent / "config" / "parameter.yaml"

def load_config(path: str = None) -> AppConfig:
    """
    YAML 파일을 읽어 Pydantic AppConfig 객체로 변환.
    Validation Error 발생 시 즉시 리포트됨.
    """
    target_path = Path(path) if path else DEFAULT_CONFIG_PATH
    
    if not target_path.exists():
        # 기본 경로 외에 현재 작업 디렉토리도 확인
        if Path("parameter.yaml").exists():
            target_path = Path("parameter.yaml")
        else:
            raise FileNotFoundError(f"Configuration file not found at: {target_path}")

    with open(target_path, "r") as f:
        data = yaml.safe_load(f)
        
    return AppConfig(**data)