from pydantic import BaseModel, Field
from typing import Dict, Any, List, Optional
# 실제 경로에 맞게 수정 필요
from pcr.components.amplicon import Amplicon

class BaseDesignInput(BaseModel):
    name: str
    target_start: int
    target_end: int
    template_sequence: str
    # [MODIFIED] process_results에서 사용되는 reference_sequence 추가
    reference_sequence: Optional[str] = None 
    template_genomic_start:  Optional[int] = None 
    template_genomic_end: Optional[int] = None 
    reference_name: str = "hg38"
    overrides: Dict[str, Any] = Field(default_factory=dict)
    
    # 파이프라인 전역 Config (PipelineConfig)
    config: Any 

class BaseDesignOutput(BaseModel):
    amplicons: List[Amplicon] = []
    status: str
    log_messages: List[str] = []
    error_msg: Optional[str] = None