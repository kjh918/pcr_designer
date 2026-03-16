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

    # 🔥 [추가] 3' 말단 윈도우 사이즈 (타겟이 프라이머 끝에서 몇 bp 이내에 와야 하는지)
    window_size_3prime: int = Field(
        default=3,
        description="Allowed window size at the 3' end to contain the target CpG (e.g.,3 means within the last 4 bases)."
    )

    # 🔥 [추가] 프라이머 서열 내부에 포함되어야 할 최소 CpG 개수
    min_cpg_count: int = Field(
        default=1,
        description="Minimum number of CpG sites required within the primer sequence."
    )

    def __init__(self, **data):
        super().__init__(**data)
        
        # 필수 키 검증
        if "M" not in self.templates or "U" not in self.templates:
            raise ValueError("MS-PCR requires both 'M' and 'U' templates.")