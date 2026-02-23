from typing import List
from ..base.qc import BaseQCExecutor
from pcr.components.amplicon import Amplicon

class QPCRQCExecutor(BaseQCExecutor):
    """
    qPCR QC 파이프라인.
    기본적으로 BaseQC(Thermo, Blast)를 수행하며, 
    어디서 탈락했는지 통계를 수집하여 Factory로 전달합니다.
    """
    
    def _setup_checkers(self):
        # 1. Base의 공통 체커(Thermo, Blast) 장착
        super()._setup_checkers()
        
        # 통계를 담을 딕셔너리 초기화 (Factory에서 읽어감)
        self.qc_stats = {}

    def execute(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        self.qc_stats["total_input"] = len(amplicons)
        
        # 1. 부모의 체인(Thermo -> Blast) 실행
        # (통과한 앰플리콘만 리턴됨)
        valid_amplicons = super().execute(amplicons)
        
        # 2. 🔥 [핵심 추가] 탈락 사유 통계 집계
        # 파이썬은 객체 참조(Reference)를 사용하므로, 
        # BaseQC가 원본 amp 객체에 남겨둔 실패 기록을 그대로 읽어올 수 있습니다.
        for amp in amplicons:
            # 통과하지 못한 객체 찾기
            if not getattr(amp, "is_qc_pass", True):
                # BaseQC의 체커들이 fail_reason을 적어두었다고 가정 (없으면 기본값)
                reason = getattr(amp, "qc_fail_reason", "BaseQC_Fail")
                self.qc_stats[reason] = self.qc_stats.get(reason, 0) + 1

        # 3. qPCR 전용 후처리 필터 (Probe 존재 여부)
        final_amplicons = []
        for amp in valid_amplicons:
            if amp.probe:
                final_amplicons.append(amp)
            else:
                # Probe가 없어서 여기서 탈락하는 경우 기록
                self.qc_stats["MissingProbe_Fail"] = self.qc_stats.get("MissingProbe_Fail", 0) + 1
                amp.is_qc_pass = False
                amp.qc_fail_reason = "Missing Probe"
                
        self.qc_stats["final_passed"] = len(final_amplicons)
        
        return final_amplicons