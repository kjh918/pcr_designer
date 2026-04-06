"""
pcr/designers/as_pcr/qc.py
AS-PCR 전용 QC 파이프라인.
세트(Set) 단위로 연대 책임을 지도록 평가하며, 
실패한 경우에도 상세 사유 확인을 위해 모든 앰플리콘을 반환합니다.
"""
from typing import List, Dict
from pcr.designers.base.qc import BaseQCExecutor
from pcr.components.amplicon import Amplicon

class ASPCRQCExecutor(BaseQCExecutor):
    """
    AS-PCR 전용 QC 실행기.
    """
    
    # 🔥 [Error Fix] 부모 클래스(BaseQCExecutor)의 호출 방식과 일치하도록 blast 인자를 추가합니다.
    def _setup_checkers(self, blast: bool = True):
        """
        AS-PCR 전용 QC 모듈 셋업.
        """
        # 부모의 물리적 특성 및 BLAST 체크 기능을 그대로 사용합니다.
        super()._setup_checkers(blast=blast)
        self.qc_stats = {}

    def execute(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        """
        1. 개별 앰플리콘 QC 수행 (Thermo, BLAST 등)
        2. 세트(Set) 단위 그룹화 및 연대 책임 판정
        """
        if not amplicons:
            return []

        self.qc_stats["total_input"] = len(amplicons)
        
        # 1. Base QC 수행 (개별 앰플리콘의 is_qc_pass, qc_log 업데이트)
        # BaseQCExecutor.execute는 내부 체커들을 돌려 리스트를 반환합니다.
        super().execute(amplicons)
        
        # 2. 세트(Set) 단위로 그룹화
        sets_dict: Dict[str, List[Amplicon]] = {}
        for amp in amplicons:
            set_id = getattr(amp, "set_id", "Unknown")
            if set_id not in sets_dict:
                sets_dict[set_id] = []
            sets_dict[set_id].append(amp)
            
            # 개별 탈락 사유 통계 수집
            if not getattr(amp, "is_qc_pass", True):
                reason = getattr(amp, "qc_log", "BaseQC_Fail")
                # 여러 사유가 섞여 있을 수 있으므로 첫 번째 주요 사유만 카운트
                primary_reason = reason.split('|')[0].strip()
                self.qc_stats[primary_reason] = self.qc_stats.get(primary_reason, 0) + 1

        # 3. 🔥 세트 연대 책임 판정 (All or Nothing)
        passed_set_count = 0
        
        for set_id, amps_in_set in sets_dict.items():
            # 판정 기준 A: 세트 내 모든 앰플리콘이 개별 QC를 통과했는가?
            is_set_pass = all(getattr(amp, "is_qc_pass", False) for amp in amps_in_set)
            
            # 판정 기준 B: AS-PCR 세트 구성(WT, ALT, WT_MM, ALT_MM) 4개가 모두 존재하는가?
            has_all_alleles = len(amps_in_set) == 4 
            
            if is_set_pass and has_all_alleles:
                passed_set_count += 1
            else:
                # 하나라도 실패하거나 누락된 경우 세트 전체를 탈락 처리
                self.qc_stats["set_failed"] = self.qc_stats.get("set_failed", 0) + 1
                for amp in amps_in_set:
                    # 이미 실패한 앰플리콘은 기존 사유 유지, 통과했던 앰플리콘만 세트 책임으로 전환
                    if getattr(amp, "is_qc_pass", True):
                        amp.is_qc_pass = False
                        fail_msg = "Set_Validation_Fail (Sibling failed QC or Missing Allele)"
                        # 기존 로그가 있다면 병합, 없으면 새로 할당
                        current_log = getattr(amp, "qc_log", "")
                        amp.qc_log = f"{current_log} | {fail_msg}".strip(" | ")

        self.qc_stats["final_passed_sets"] = passed_set_count
        
        # 결과 리포팅을 위해 모든 앰플리콘 리스트를 그대로 반환 (is_qc_pass 플래그만 수정됨)
        return amplicons