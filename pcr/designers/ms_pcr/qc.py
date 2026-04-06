"""
pcr/designers/ms_pcr/qc.py
MS-PCR 전용 QC 파이프라인.
M-Allele과 U-Allele 앰플리콘이 한 쌍(Set)으로 묶여 연대 책임을 지도록 평가합니다.
"""
from typing import List, Dict
from pcr.designers.base.qc import BaseQCExecutor
from pcr.components.amplicon import Amplicon

class MSPCRQCExecutor(BaseQCExecutor):
    """
    MS-PCR 전용 QC 파이프라인.
    """
    
    def _setup_checkers(self, blast: bool = False):
        # 🔥 Bisulfite 변환 프라이머는 일반 게놈(hg38) BLAST 시 100% 탈락하므로,
        # MS-PCR 전용 변환 DB가 구축되기 전까지는 기본적으로 BLAST 체커를 끕니다.
        super()._setup_checkers(blast=False) 
        self.qc_stats = {}

    def execute(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        self.qc_stats["total_input"] = len(amplicons)
        
        # 1. Base QC 수행 (Tm, GC, Dimer 등 열역학 검사)
        super().execute(amplicons)
        
        # 2. 세트(Set) 단위로 그룹화 (M과 U)
        sets_dict: Dict[str, List[Amplicon]] = {}
        for amp in amplicons:
            set_id = getattr(amp, "set_id", "Unknown")
            if set_id not in sets_dict:
                sets_dict[set_id] = []
            sets_dict[set_id].append(amp)
            
            # 개별 탈락 사유 통계 수집
            if not getattr(amp, "is_qc_pass", True):
                reason = getattr(amp, "qc_log", "BaseQC_Fail")
                self.qc_stats[reason] = self.qc_stats.get(reason, 0) + 1

        # 3. 🔥 세트 연대 책임 판정 (All or Nothing)
        passed_set_count = 0
        
        for set_id, amps_in_set in sets_dict.items():
            # MS-PCR은 한 세트에 반드시 "M"과 "U" 2개가 존재해야 함
            is_set_pass = all(getattr(amp, "is_qc_pass", False) for amp in amps_in_set)
            
            # Allele 검증: M과 U가 모두 있는지 확인
            allele_types = {getattr(amp, "allele_type", "unknown") for amp in amps_in_set}
            has_all_alleles = ("M" in allele_types and "U" in allele_types)
            
            if is_set_pass and has_all_alleles:
                passed_set_count += 1
            else:
                self.qc_stats["set_failed"] = self.qc_stats.get("set_failed", 0) + 1
                # 🚨 한쪽이라도 탈락했다면, 통과한 나머지 형제 앰플리콘도 강등(False)시킴
                for amp in amps_in_set:
                    if getattr(amp, "is_qc_pass", True):
                        amp.is_qc_pass = False
                        amp.qc_log = "Set_Validation_Fail (Sibling M or U allele failed QC)"

        self.qc_stats["final_passed_sets"] = passed_set_count
        
        # 🔥 실패 사유 리포팅을 위해 Pydantic 드롭 없이 원본 리스트 전체 반환
        return amplicons