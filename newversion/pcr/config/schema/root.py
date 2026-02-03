from typing import Dict, Optional, Any, List
from pydantic import BaseModel, Field

# 하위 스키마 임포트
from .references import ReferenceConfig
from .pcr import PCRParams
#from .pcr.components import Amplicon
from ...components.amplicon import Amplicon
from .qc import QCCriteria, QCToolsConfig

# ----------------------------------------------------------------
# 1. From system.yaml (환경 설정)
# ----------------------------------------------------------------

class AppConfig(BaseModel):
    # (1) References: 종(Species) 이름이 키가 됨 (hg19, hg38...)
    references: Dict[str, ReferenceConfig]    
    # (2) QC Tools: 실행 파일 경로들
    qc_tools: QCToolsConfig
    # (3) PCR Parameters: Primer/Probe 디자인 조건
    pcr_params: PCRParams
    # (4) QC Criteria: 통과/탈락 기준
    qc_criteria: QCCriteria

# =============================================================================
# 2. 실행 입력 (Runtime Input) - 서열, 좌표 등
# =============================================================================
class BaseDesignInput(BaseModel):
    """Designer에게 일을 시킬 때 전달하는 봉투"""
    # 필수 데이터
    name: str = Field(..., description="Unique name for this primer design task")
    template_sequence: str = Field(..., description="Target Template Sequence")
    target_start: int = Field(..., ge=0, description="0-based Relative Start")
    target_end: int = Field(..., ge=0, description="0-based Relative End")
    
    # 설정 정보 주입
    config: AppConfig
    
    # 메타데이터 (참조 유전체 정보 등)
    reference_name: str = "hg38" # hg38, hg19 등 정보만 유지
    reference_sequence: Optional[str] = None
    
    # 런타임 오버라이드
    overrides: Dict[str, Any] = Field(default_factory=dict)

    #@root_validator(pre=True)
    def check_range(cls, values):
        s, e = values.get('target_start'), values.get('target_end')
        if s is not None and e is not None and e <= s:
             raise ValueError(f"Target end ({e}) must be > start ({s})")
        return values
		
# =============================================================================
# 3. 실행 결과 (Runtime Output) - 결과 앰플리콘 리스트 등
# =============================================================================
class BaseDesignOutput(BaseModel):
    """Designer가 업무를 마치고 반환하는 결과 보고서"""
    amplicons: List[Amplicon] = Field(default_factory=list)
    status: str = Field(..., description="success, no_candidates, error")
    error_msg: Optional[str] = None
    log_messages: List[str] = Field(default_factory=list)

    class Config:
        arbitrary_types_allowed = True

    @property
    def passed_count(self) -> int:
        return len([a for a in self.amplicons if a.is_qc_pass])