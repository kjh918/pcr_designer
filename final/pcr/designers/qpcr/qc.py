"""
pcr/designers/qpcr/qc.py
TaqMan qPCR 전용 QC 파이프라인.
Base QC(Thermo, Blast)를 통과한 앰플리콘들을 대상으로 qPCR 특화 검증을 추가 수행합니다.
"""
from typing import List
from pcr.designers.base.qc import BaseQCExecutor
from pcr.components.amplicon import Amplicon

# =================================================================
# [모듈화된 체커(Checker) 클래스들]
# 스키마의 하위 그룹(primer, probe 등)을 주입받아 독립적으로 검증을 수행합니다.
# 다른 PCR 기법에서도 원한다면 이 클래스들을 가져다 쓸 수 있습니다.
# =================================================================

class QPCRPrimerChecker:
    """Primer 세트의 물리적 조건(ΔTm 등)을 검사하는 재사용 가능한 모듈"""
    def __init__(self, primer_criteria):
        # 스키마의 'primer' 관련 조건 세트만 주입받음
        self.criteria = primer_criteria
        
    def run(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        max_tm_diff = getattr(self.criteria, "max_diff_tm", 3.0)
        
        for amp in amplicons:
            # 이미 이전 체커에서 심각한 결격 사유로 탈락한 경우 스킵 (선택사항)
            # 여기서는 UI 표시를 위해 실패했더라도 사유를 계속 누적시킵니다.
            tm_diff = abs(amp.forward.tm - amp.reverse.tm)
            
            if tm_diff > max_tm_diff:
                amp.is_qc_pass = False
                amp.qc_log = getattr(amp, "qc_log", "") + f" [Primer ΔTm {round(tm_diff,1)} > {max_tm_diff}]"
                
        return amplicons

class QPCRProbeChecker:
    """TaqMan Probe의 품질 조건(5' G, Poly-G, Tm 등)을 검사하는 재사용 가능한 모듈"""
    def __init__(self, probe_criteria):
        # 스키마의 'probe' 관련 조건 세트만 주입받음
        self.criteria = probe_criteria
        
    def run(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        avoid_5g = getattr(self.criteria, "avoid_5_prime_g", True)
        max_poly_g = getattr(self.criteria, "max_probe_poly_g", 3)
        
        for amp in amplicons:
            if not amp.probe:
                continue
                
            p_seq = amp.probe.sequence.upper()
            qc_log = getattr(amp, "qc_log", "")
            is_pass = getattr(amp, "is_qc_pass", True)
            
            # 규칙 A: 5' 말단 G 금지 (형광 염료 Quenching 방지)
            if avoid_5g and p_seq.startswith('G'):
                is_pass = False
                qc_log += " [Probe 5' starts with 'G' (Quenching risk)]"
            
            # 규칙 B: Probe Tm은 두 Primer의 Tm 중 높은 것보다 높아야 함
            max_pr_tm = max(amp.forward.tm, amp.reverse.tm)
            if amp.probe.tm <= max_pr_tm:
                is_pass = False
                qc_log += f" [Probe Tm({round(amp.probe.tm,1)}) ≤ Max Primer Tm({round(max_pr_tm,1)})]"
            
            # 규칙 C: Poly-G 검사 (G가 연속 4개 이상 나오는 것 방지 등)
            if "G" * (max_poly_g + 1) in p_seq:
                is_pass = False
                qc_log += f" [Probe contains Poly-G (>{max_poly_g})]"
            
            # 상태 업데이트
            amp.is_qc_pass = is_pass
            amp.qc_log = qc_log
                
        return amplicons

# =================================================================
# [메인 Executor]
# =================================================================

class QPCRQCExecutor(BaseQCExecutor):
    """
    qPCR 전용 QC 파이프라인.
    """
    
    def _setup_checkers(self):
        # 1. Base 체커 세팅 (Thermo, Blast)
        use_blast = getattr(self.config.system, "blast_db_path", None) is not None
        super()._setup_checkers(blast=use_blast)
        
        # 2. qPCR 전용 체커 추가
        # 스키마(config.qc_criteria)에서 각각 필요한 그룹 세트만 잘라서 주입합니다.
        qc_criteria = self.config.qc_criteria
        
        if hasattr(qc_criteria, "primer"):
            self.checkers.append(QPCRPrimerChecker(qc_criteria.primer))
            
        if hasattr(qc_criteria, "probe"):
            self.checkers.append(QPCRProbeChecker(qc_criteria.probe))

    def execute(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        # 부모의 execute()를 호출하면, self.checkers에 담긴 모든 체커들이 
        # 순서대로 앰플리콘을 검증(run)하고 통과시켜줍니다.
        evaluated_amps = super().execute(amplicons)
        
        # UI 표기를 위해 공백 정리만 수행
        for amp in evaluated_amps:
            amp.qc_log = getattr(amp, "qc_log", "").strip()
            
        return evaluated_amps