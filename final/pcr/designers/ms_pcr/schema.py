"""
pcr/designers/ms_pcr/schema.py
MS-PCR 전용 입력 스키마 정의.
"""
from typing import Dict, List, Optional
from pydantic import Field

from ..base.schema import BaseDesignInput

class MSPCRDesignInput(BaseDesignInput):
    """
    MS-PCR 설계를 위한 입력 데이터 모델.
    BaseDesignInput의 기본 속성(name, target_start 등)을 상속받고,
    Bisulfite Conversion 관련 데이터를 추가로 받습니다.
    """
    
    # M-Allele (메틸화), U-Allele (비메틸화) 서열을 담은 딕셔너리
    templates: Dict[str, str] = Field(
        ..., 
        description="Dictionary containing 'M' and 'U' converted sequences"
    )
    
    # 타겟으로 삼은 CpG의 1-based 인덱스 리스트
    target_cpg_indices: List[int] = Field(
        default_factory=list,
        description="List of 1-based indices of target CpG sites in the raw sequence"
    )

    def __init__(self, **data):
        super().__init__(**data)
        
        # 필수 키 검증
        if "M" not in self.templates or "U" not in self.templates:
            raise ValueError("MS-PCR requires both 'M' and 'U' templates.")