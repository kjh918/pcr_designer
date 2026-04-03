"""
pcr/designers/qpcr/qc.py
TaqMan qPCR 전용 QC 파이프라인.
프라이머와 프로브의 개별 물리적 특성을 검증하고, 
BLAST 결과를 바탕으로 오프타겟 시그널을 분리하여 리포팅합니다.
"""
import copy
from typing import List, Dict, Any
from pcr.designers.base.qc import BaseQCExecutor
from pcr.components.amplicon import Amplicon
from pcr.config.schema.qc import AmpliconQCStatus

class QPCRPrimerChecker:
    def __init__(self, qc_criteria):
        # 전체 qc_criteria를 받아와서 dG 한계값(limit)을 읽을 수 있도록 설정
        self.criteria = qc_criteria
        self.primer_crit = getattr(qc_criteria, "primer", qc_criteria)
        
    def run(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        max_tm_diff = getattr(self.primer_crit, "max_diff_tm", 3.0)
        
        # 2차 구조(dG) 한계치 가져오기
        hp_limit = getattr(self.criteria, "hairpin_min_dg", -5.0)
        hd_limit = getattr(self.criteria, "homodimer_min_dg", -6.0)
        he_limit = getattr(self.criteria, "heterodimer_min_dg", -6.0)

        for amp in amplicons:
            if not hasattr(amp, "qc_status") or amp.qc_status is None:
                amp.qc_status = AmpliconQCStatus()
            
            is_pass = True
            msgs = []

            # 1. Primer ΔTm 검사
            tm_diff = abs(amp.forward.tm - amp.reverse.tm)
            if tm_diff > max_tm_diff:
                is_pass = False
                msgs.append(f"Primer ΔTm ({round(tm_diff,1)}℃) > {max_tm_diff}℃")

            # 2. 2차 구조 (Secondary Structure dG) 검사
            f_hp = getattr(amp.forward, "hairpin_dg", 0.0)
            r_hp = getattr(amp.reverse, "hairpin_dg", 0.0)
            f_hd = getattr(amp.forward, "homodimer_dg", 0.0)
            r_hd = getattr(amp.reverse, "homodimer_dg", 0.0)
            
            # Thermo 모듈에서 계산된 Heterodimer 데이터 추출
            t_metrics = {}
            if hasattr(amp.qc_status, "modules") and "thermo" in amp.qc_status.modules:
                t_metrics = amp.qc_status.modules["thermo"].metrics
            fr_he = t_metrics.get("hetero_fr_dg", 0.0)

            if f_hp < hp_limit:
                is_pass = False; msgs.append(f"Fwd Hairpin ({round(f_hp,2)}) < {hp_limit}")
            if r_hp < hp_limit:
                is_pass = False; msgs.append(f"Rev Hairpin ({round(r_hp,2)}) < {hp_limit}")
            if f_hd < hd_limit:
                is_pass = False; msers.append(f"Fwd Homodimer ({round(f_hd,2)}) < {hd_limit}")
            if r_hd < hd_limit:
                is_pass = False; msgs.append(f"Rev Homodimer ({round(r_hd,2)}) < {hd_limit}")
            if fr_he < he_limit:
                is_pass = False; msgs.append(f"F/R Heterodimer ({round(fr_he,2)}) < {he_limit}")

            amp.qc_status.add_result(
                module_name="qpcr_primer", 
                is_pass=is_pass, 
                messages=msgs, 
                metrics={"tm_diff": round(tm_diff, 2)}
            )
            amp.is_qc_pass = amp.qc_status.is_pass
            amp.qc_log = " | ".join(amp.qc_status.fail_reasons)
            
        return amplicons


class QPCRProbeChecker:
    def __init__(self, qc_criteria):
        self.criteria = qc_criteria
        self.probe_crit = getattr(qc_criteria, "probe", qc_criteria)
        
    def run(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        avoid_5g = getattr(self.probe_crit, "avoid_5_prime_g", True)
        max_poly_g = getattr(self.probe_crit, "max_probe_poly_g", 3)
        min_probe_diff = getattr(self.probe_crit, "min_primer_probe_tm_diff", 5.0)
        max_probe_diff = getattr(self.probe_crit, "max_primer_probe_tm_diff", 10.0)
        
        hp_limit = getattr(self.criteria, "hairpin_min_dg", -5.0)
        hd_limit = getattr(self.criteria, "homodimer_min_dg", -6.0)
        he_limit = getattr(self.criteria, "heterodimer_min_dg", -6.0)

        for amp in amplicons:
            if not amp.probe: continue
            if not hasattr(amp, "qc_status") or amp.qc_status is None: 
                amp.qc_status = AmpliconQCStatus()
                
            p_seq = amp.probe.sequence.upper()
            is_pass = True
            msgs = []
            
            # 1. 기본 Probe 규칙 검사
            if avoid_5g and p_seq.startswith('G'):
                is_pass = False; msgs.append("Probe 5' starts with 'G'")
            
            max_pr_tm = max(amp.forward.tm, amp.reverse.tm)
            target_min_tm = max_pr_tm + min_probe_diff
            target_max_tm = max_pr_tm + max_probe_diff
            
            if not (target_min_tm <= amp.probe.tm <= target_max_tm):
                is_pass = False
                msgs.append(f"Probe Tm ({round(amp.probe.tm,1)}℃) out of range. Expected: {round(target_min_tm,1)}~{round(target_max_tm,1)}℃")
            
            if "G" * (max_poly_g + 1) in p_seq:
                is_pass = False; msgs.append(f"Probe contains Poly-G (>{max_poly_g})")
                
            # 2. Probe 2차 구조 (Secondary Structure dG) 검사
            p_hp = getattr(amp.probe, "hairpin_dg", 0.0)
            p_hd = getattr(amp.probe, "homodimer_dg", 0.0)
            
            t_metrics = {}
            if hasattr(amp.qc_status, "modules") and "thermo" in amp.qc_status.modules:
                t_metrics = amp.qc_status.modules["thermo"].metrics
                
            fp_he = t_metrics.get("hetero_fp_dg", 0.0)
            rp_he = t_metrics.get("hetero_rp_dg", 0.0)

            if p_hp < hp_limit:
                is_pass = False; msgs.append(f"Probe Hairpin ({round(p_hp,2)}) < {hp_limit}")
            if p_hd < hd_limit:
                is_pass = False; msgs.append(f"Probe Homodimer ({round(p_hd,2)}) < {hd_limit}")
            if fp_he < he_limit:
                is_pass = False; msgs.append(f"F/P Heterodimer ({round(fp_he,2)}) < {he_limit}")
            if rp_he < he_limit:
                is_pass = False; msgs.append(f"R/P Heterodimer ({round(rp_he,2)}) < {he_limit}")

            amp.qc_status.add_result(
                module_name="qpcr_probe", 
                is_pass=is_pass, 
                messages=msgs, 
                metrics={"avoid_5g": avoid_5g, "max_poly_g": max_poly_g}
            )
            amp.is_qc_pass = amp.qc_status.is_pass
            amp.qc_log = " | ".join(amp.qc_status.fail_reasons)
            
        return amplicons

class QPCRQCExecutor(BaseQCExecutor):
    """
    qPCR 전용 QC 실행기.
    """
    def _setup_checkers(self, blast: bool = True):
        """
        qPCR 전용 QC 모듈을 셋업합니다.
        BaseQCExecutor에서 결정된 blast 실행 여부를 그대로 수용합니다.
        """
        # 부모 클래스의 기본 체커(Thermo, Blast) 등록
        super()._setup_checkers(blast=blast)
        
        qc_criteria = self.config.qc_criteria
        # qPCR 전용 정밀 체커 추가
        self.checkers.append(QPCRPrimerChecker(qc_criteria))
        self.checkers.append(QPCRProbeChecker(qc_criteria))

    def execute(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        """
        1. BaseQC(Thermo, Blast) + qPCR 전용 체커 실행
        2. BLAST 결과에 따른 앰플리콘 분리 (Explode)
        """
        # 모든 체커 순차 실행
        evaluated_amps = super().execute(amplicons)
        
        expanded_amplicons = []
        for amp in evaluated_amps:
            blast_details = getattr(amp, "blast_stats", {})
            
            # BLAST 결과 시그널이 하나라도 있는 경우 분리 로직 수행
            if blast_details and blast_details.get("total_signal_count", 0) > 0:
                rank = 1
                aln_base = getattr(amp, "alignment_visual", [])
                
                # A. Target 시그널 분리
                for sig in blast_details.get("target_signals", []):
                    cloned = copy.deepcopy(amp)
                    cloned.id = f"{amp.id}_Target_{rank}"
                    cloned.product_size = sig.get("product_size", amp.product_size)
                    cloned.genomic_pos = sig.get("location", "Unknown")
                    # 원본 정렬 정보 위에 BLAST 상세 정렬 정보 추가
                    cloned.alignment_visual = aln_base + [""] + [sig.get("unified_text_block", "")]
                    expanded_amplicons.append(cloned)
                    rank += 1
                    
                # B. Off-Target 시그널 분리 (FAIL 처리)
                for sig in blast_details.get("off_target_signals", []):
                    cloned = copy.deepcopy(amp)
                    cloned.id = f"{amp.id}_OffTarget_{rank}"
                    cloned.product_size = sig.get("product_size", "N/A")
                    cloned.genomic_pos = sig.get("location", "Unknown")
                    cloned.alignment_visual = [sig.get("unified_text_block", "")]
                    cloned.is_qc_pass = False
                    cloned.qc_log = f"BLAST: Off-Target (Probe binds here!) | {cloned.qc_log}".strip(" |")
                    expanded_amplicons.append(cloned)
                    rank += 1
                    
                # C. Noise (앰플리콘만 증폭) 시그널 분리 (FAIL 처리)
                for sig in blast_details.get("amplification_only", []):
                    cloned = copy.deepcopy(amp)
                    cloned.id = f"{amp.id}_Noise_{rank}"
                    cloned.product_size = sig.get("product_size", "N/A")
                    cloned.genomic_pos = sig.get("location", "Unknown")
                    cloned.alignment_visual = [sig.get("unified_text_block", "")]
                    cloned.is_qc_pass = False
                    cloned.qc_log = f"BLAST: Amplification Only | {cloned.qc_log}".strip(" |")
                    expanded_amplicons.append(cloned)
                    rank += 1
            else:
                # BLAST를 돌렸으나 결과가 전혀 없는 경우 혹은 BLAST를 수행하지 않은 경우
                pcr_params = getattr(self.config, "pcr_params", None)
                ref_name = getattr(pcr_params, "reference_name", "none") if pcr_params else "none"
                
                # 유전체가 hg38 등으로 지정되었으나 BLAST 시그널이 없는 경우에만 실패 처리
                if str(ref_name).lower() != "none" and blast_details:
                    amp.is_qc_pass = False
                    amp.qc_log = f"BLAST: No target found in reference genome. | {amp.qc_log}".strip(" |")
                
                expanded_amplicons.append(amp)

        return expanded_amplicons