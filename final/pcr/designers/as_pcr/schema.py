from typing import Dict, Optional, List
from pydantic import Field
from pcr.designers.base.schema import BaseDesignInput, BaseDesignOutput
from pcr.components.amplicon import Amplicon

class ASPCRAmplicon(Amplicon):
    """
    AS-PCR 전용 앰플리콘 모델. 
    일반 Amplicon에 AS-PCR 분석 및 시각화용 메타데이터 필드를 추가합니다.
    """
    allele_type: Optional[str] = None
    set_id: Optional[str] = None
    fixed_prime: Optional[str] = None
    alignment_visual: Optional[List[str]] = None

class ASPCRDesignInput(BaseDesignInput):
    """
    AS-PCR 전용 입력 스키마.
    BaseDesignInput을 상속받으며, 4가지 템플릿 딕셔너리와 앵커 방향이 추가됩니다.
    """
    templates: Dict[str, str] = Field(
        ..., 
        description="WT, ALT, WT_MM, ALT_MM 서열을 포함하는 딕셔너리"
    )
    fixed_prime: str = Field(
        default="forward", 
        description="SNP 위치에 3' 말단을 고정할 프라이머 방향 ('forward' 또는 'reverse')"
    )

class ASPCRDesignOutput(BaseDesignOutput):
    """
    AS-PCR 전용 출력 스키마.
    향후 세트(Set) 단위의 그룹화 정보나 AS-PCR 특화 메타데이터를 담을 때 확장합니다.
    """
    pass