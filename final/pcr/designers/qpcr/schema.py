"""
pcr/designers/qpcr/schema.py
TaqMan qPCR 전용 코어 엔진 데이터 검증 스키마.
"""
from typing import Optional
from pydantic import Field
from pcr.designers.base.schema import BaseDesignInput, BaseDesignOutput

class QPCRDesignInput(BaseDesignInput):
    """
    qPCR 설계를 위한 코어 엔진 입력 데이터 모델.
    부모 클래스(BaseDesignInput)로부터 대괄호 파싱(template_sequence, target_start, target_end) 
    기능을 완벽하게 상속받으며, 아래에 qPCR만의 고유한 특성을 추가합니다.
    """
    
    # 🔥 qPCR 전용 파라미터 확장
    require_probe: bool = Field(
        default=True,
        description="TaqMan qPCR의 경우 내부 프로브(Internal Probe) 설계가 필수인지 여부"
    )
    
    
class QPCRDesignOutput(BaseDesignOutput):
    """
    qPCR 설계 결과 모델.
    부모(BaseDesignOutput)의 amplicons, status, log_messages 속성을 상속받습니다.
    """
    
    # 🔥 qPCR 전용 결과 메타데이터 확장
    probe_designed_count: int = Field(
        default=0, 
        description="설계에 성공하여 프로브가 포함된 최종 앰플리콘 세트의 개수"
    )