"""
pcr/designers/qpcr/qc.py
TaqMan qPCR 전용 QC 파이프라인
"""
from typing import List
from ..base.qc import BaseQCExecutor
from pcr.components.amplicon import Amplicon

class QPCRQCExecutor(BaseQCExecutor):
    """
    [MODIFIED] qPCR QC 파이프라인.
    기본적으로 BaseQC(Thermo, Blast)를 수행하며, 
    ThermoChecker 내부에서 Probe와의 Heterodimer 간섭을 자동으로 잡아냅니다.
    """
    
    def _setup_checkers(self):
        # 1. Base의 공통 체커(Thermo, Blast) 장착
        super()._setup_checkers()
        
        # 2. [추가 확장 포인트] qPCR 전용 커스텀 체커가 있다면 여기에 추가
        # 예: self.checkers.append(QPCRSpecificDyeChecker())

    def execute(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        # 부모의 체인(Thermo -> Blast) 실행
        valid_amplicons = super().execute(amplicons)
        
        # [추가 확장 포인트] 체커(Class)로 빼기 애매한 간단한 후처리 필터가 있다면 여기서 수행
        final_amplicons = []
        for amp in valid_amplicons:
            # 예: qPCR은 Probe가 무조건 있어야만 최종 Pass 처리
            if amp.probe:
                final_amplicons.append(amp)
                
        return final_amplicons