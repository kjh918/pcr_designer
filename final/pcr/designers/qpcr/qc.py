"""
pcr/designers/qpcr/qc.py
TaqMan qPCR 전용 QC 파이프라인.
"""
import copy
from typing import List
from pcr.designers.base.qc import BaseQCExecutor
from pcr.components.amplicon import Amplicon
from pcr.config.schema.qc import AmpliconQCStatus


import copy
from typing import List
from pcr.designers.base.qc import BaseQCExecutor
from pcr.components.amplicon import Amplicon
from pcr.config.schema.qc import AmpliconQCStatus

class QPCRPrimerChecker:
    def __init__(self, qc_criteria):
        # 💡 전체 qc_criteria를 받아와서 dG 한계값(limit)을 읽을 수 있도록 수정
        self.criteria = qc_criteria
        self.primer_crit = getattr(qc_criteria, "primer", qc_criteria)
        
    def run(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        max_tm_diff = getattr(self.primer_crit, "max_diff_tm", 3.0)
        
        # 2차 구조(dG) 한계치 가져오기 (없으면 기본값 적용)
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
            # (dG는 음수이므로, 계산값이 한계값보다 작으면 더 강하게 결합한다는 뜻 = FAIL)
            f_hp = getattr(amp.forward, "hairpin_dg", 0.0)
            r_hp = getattr(amp.reverse, "hairpin_dg", 0.0)
            f_hd = getattr(amp.forward, "homodimer_dg", 0.0)
            r_hd = getattr(amp.reverse, "homodimer_dg", 0.0)
            
            # Heterodimer는 두 프라이머 간의 상호작용이므로 thermo metrics에서 추출
            t_metrics = {}
            if hasattr(amp.qc_status, "modules") and "thermo" in amp.qc_status.modules:
                t_metrics = amp.qc_status.modules["thermo"].metrics
            fr_he = t_metrics.get("hetero_fr_dg", 0.0)

            if f_hp < hp_limit:
                is_pass = False; msgs.append(f"Fwd Hairpin ({round(f_hp,2)}) < {hp_limit}")
            if r_hp < hp_limit:
                is_pass = False; msgs.append(f"Rev Hairpin ({round(r_hp,2)}) < {hp_limit}")
            if f_hd < hd_limit:
                is_pass = False; msgs.append(f"Fwd Homodimer ({round(f_hd,2)}) < {hd_limit}")
            if r_hd < hd_limit:
                is_pass = False; msgs.append(f"Rev Homodimer ({round(r_hd,2)}) < {hd_limit}")
            if fr_he < he_limit:
                is_pass = False; msgs.append(f"F/R Heterodimer ({round(fr_he,2)}) < {he_limit}")

            amp.qc_status.add_result(module_name="qpcr_primer", is_pass=is_pass, messages=msgs, metrics={"tm_diff": round(tm_diff, 2)})
            amp.is_qc_pass = amp.qc_status.is_pass
            amp.qc_log = " | ".join(amp.qc_status.fail_reasons)
            
        return amplicons


class QPCRProbeChecker:
    def __init__(self, qc_criteria):
        # 💡 전체 qc_criteria를 받아와서 dG 한계값(limit)을 읽을 수 있도록 수정
        self.criteria = qc_criteria
        self.probe_crit = getattr(qc_criteria, "probe", qc_criteria)
        
    def run(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        avoid_5g = getattr(self.probe_crit, "avoid_5_prime_g", True)
        max_poly_g = getattr(self.probe_crit, "max_probe_poly_g", 3)
        
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
            if amp.probe.tm <= max_pr_tm:
                is_pass = False; msgs.append(f"Probe Tm({round(amp.probe.tm,1)}) ≤ Max Primer Tm({round(max_pr_tm,1)})")
            
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

            amp.qc_status.add_result(module_name="qpcr_probe", is_pass=is_pass, messages=msgs, metrics={"avoid_5g": avoid_5g, "max_poly_g": max_poly_g})
            amp.is_qc_pass = amp.qc_status.is_pass
            amp.qc_log = " | ".join(amp.qc_status.fail_reasons)
            
        return amplicons


class QPCRQCExecutor(BaseQCExecutor):
    def _setup_checkers(self):
        # BaseQCExecutor가 알아서 BlastSpecificityChecker를 추가하도록 설정
        use_blast = False
        if hasattr(self.config, "system") and hasattr(self.config.system, "paths"):
            use_blast = getattr(self.config.system.paths, "blast_db_path", None) is not None
        super()._setup_checkers(blast=use_blast)
        
        qc_criteria = self.config.qc_criteria
        if hasattr(qc_criteria, "primer"): self.checkers.append(QPCRPrimerChecker(qc_criteria.primer))
        if hasattr(qc_criteria, "probe"): self.checkers.append(QPCRProbeChecker(qc_criteria.probe))

    def execute(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        # 1. BaseQC + BLAST + qPCR 체커 모두 실행
        evaluated_amps = super().execute(amplicons)
        
        # 2. BLAST 결과에 따라 앰플리콘을 여러 개로 증식(Explode)시킴
        expanded_amplicons = []
        for amp in evaluated_amps:
            amp.qc_log = getattr(amp, "qc_log", "").strip()
            blast_details = getattr(amp, "blast_stats", {})
            
            if blast_details and blast_details.get("total_signal_count", 0) > 0:
                rank = 1
                aln_base = getattr(amp, "alignment_visual", [])
                
                for sig in blast_details.get("target_signals", []):
                    cloned = copy.deepcopy(amp)
                    cloned.id = f"{amp.id}_Target_{rank}"
                    cloned.product_size = sig.get("product_size", "N/A")
                    cloned.genomic_pos = sig.get("location", "Unknown")
                    cloned.alignment_visual = aln_base + [""] + [sig.get("unified_text_block", "")]
                    expanded_amplicons.append(cloned)
                    rank += 1
                    
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
                # BLAST는 돌렸는데 타겟이 없는 경우
                use_blast = getattr(self.config.system.paths, "blast_db_path", None) is not None if hasattr(self.config.system, "paths") else False
                if use_blast:
                    amp.is_qc_pass = False
                    amp.qc_log = f"BLAST: No target found in reference genome. | {amp.qc_log}".strip(" |")
                expanded_amplicons.append(amp)

        return expanded_amplicons