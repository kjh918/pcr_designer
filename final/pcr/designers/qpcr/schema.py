"""
pcr/designers/qpcr/schema.py
TaqMan qPCR 전용 데이터 검증 스키마.
"""
from pcr.designers.base.schema import BaseDesignInput, BaseDesignOutput
from typing import Dict, Optional, List
from pydantic import Field
from pcr.designers.base.schema import BaseDesignInput, BaseDesignOutput
from pcr.components.amplicon import Amplicon # 상속을 위해 임포트

# ... (기존 ASPCRDesignInput 코드는 그대로 유지) ...

class ASPCRAmplicon(Amplicon):
    """
    AS-PCR 전용 앰플리콘 모델. 
    일반 Amplicon에 AS-PCR 분석 및 시각화용 메타데이터 필드를 추가합니다.
    """
    allele_type: Optional[str] = None
    set_id: Optional[str] = None
    fixed_prime: Optional[str] = None
    alignment_visual: Optional[List[str]] = None
    
class QPCRDesignInput(BaseDesignInput):
    """
    [MODIFIED] 확장 포인트 마련:
    현재는 BaseDesignInput의 구조를 그대로 사용하지만, 
    추후 qPCR에만 필요한 입력값(예: fluorophore 선호도 등)이 생기면 여기에 추가합니다.
    """
    pass

class QPCRDesignOutput(BaseDesignOutput):
    """
    [MODIFIED] 확장 포인트 마련:
    마찬가지로 qPCR 전용 결과(예: 특정 Probe의 형광 효율 계산값 등)를 
    내보내야 할 경우 이 클래스를 확장하여 사용합니다.
    """
    pass