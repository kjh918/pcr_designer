from pydantic import BaseModel
from typing import Dict
from .pcr import PCRParams
from .qc import QCCriteria
from .references import SystemConfig, ReferenceConfig

class PipelineConfig(BaseModel):
    """
    [MODIFIED] system과 references 스키마가 추가되어 전체 파이프라인의 단일 진실 공급원이 됩니다.
    """
    system: SystemConfig
    references: Dict[str, ReferenceConfig]
    pcr_params: PCRParams
    qc_criteria: QCCriteria

    def get_reference(self, ref_name: str) -> ReferenceConfig:
        """안전하게 특정 레퍼런스(hg38 등)의 경로를 반환하는 헬퍼 메서드"""
        if ref_name not in self.references:
            raise ValueError(f"Reference '{ref_name}' is not defined in system.yaml.")
        return self.references[ref_name]