from typing import List, Dict
from pcr.designers.base.qc import BaseQCExecutor
from pcr.components.amplicon import Amplicon

class ASPCRQCExecutor(BaseQCExecutor):
    """
    AS-PCR 전용 QC 파이프라인.
    세트(Set) 단위로 연대 책임을 지도록 평가합니다.
    """
    
    def _setup_checkers(self):
        super()._setup_checkers(blast=True)
        self.qc_stats = {}

    def execute(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        self.qc_stats["total_input"] = len(amplicons)
        
        # 1. 개별 앰플리콘에 대해 BaseQC(Thermo, Blast 등) 수행
        # (주의: super().execute()는 통과한 객체만 리턴하지만, 
        # 원본 amplicons 리스트 내 객체들의 is_qc_pass 상태는 업데이트되어 있음)
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
                reason = getattr(amp, "qc_fail_reason", "BaseQC_Fail")
                self.qc_stats[reason] = self.qc_stats.get(reason, 0) + 1

        # 3. 🔥 세트 연대 책임 판정 (All or Nothing)
        final_amplicons = []
        for set_id, amps_in_set in sets_dict.items():
            # 4가지(WT, ALT, WT_MM, ALT_MM)가 모두 존재하고, 모두 Pass 했는가?
            is_set_pass = all(getattr(amp, "is_qc_pass", False) for amp in amps_in_set)
            has_all_alleles = len(amps_in_set) == 4 
            
            if is_set_pass and has_all_alleles:
                # 완벽한 세트만 최종 통과 명단에 올림
                final_amplicons.extend(amps_in_set)
                self.qc_stats["set_passed"] = self.qc_stats.get("set_passed", 0) + 1
            else:
                self.qc_stats["set_failed"] = self.qc_stats.get("set_failed", 0) + 1
                # 통과했던 형제 앰플리콘들도 세트 탈락 사유를 적어 강등시킴
                for amp in amps_in_set:
                    if getattr(amp, "is_qc_pass", True):
                        amp.is_qc_pass = False
                        amp.qc_fail_reason = "Set_Validation_Fail (Sibling failed QC or Missing Allele)"

        self.qc_stats["final_passed_amplicons"] = len(final_amplicons)
        
        # 💡 옵션: 프론트엔드에 "실패한 세트"도 보여주고 싶다면 
        # return amplicons 를 하셔도 됩니다. (보통은 final_amplicons 반환)
        return final_amplicons