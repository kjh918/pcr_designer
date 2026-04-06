import pandas as pd
from typing import List, Optional, Dict, Any
from Bio.Seq import Seq
# 1. Config & Schema Import
from ...config.schema.app import AppConfig
from ...config.schema.qc import AmpliconQCStatus

# 2. Components
from ...components.amplicon import Amplicon

# 3. QC Modules
from .thermo import ThermoChecker
from .blast import BlastSpecificityChecker


class QCExecutor:
    """
    QC 실행기: Thermo -> Specificity (BLAST Only)
    """
    def __init__(self, config: AppConfig):
        self.config = config
        
        # 1. Thermo Checker (물리적 성질)
        self.thermo_checker = ThermoChecker(config.qc_criteria)
        
        # 2. BLAST Specificity Checker (특이성 검사)
        # 이제 이 모듈이 Off-target 판별의 모든 책임을 집니다.
        self.spec_checker = BlastSpecificityChecker(config)
        

    def run_qc(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        processed_amplicons = []

        for amp in amplicons:
            # -----------------------------------------------------
            # Step 1: Thermo QC
            # -----------------------------------------------------
            thermo_res = self.thermo_checker.check(amp)
            
            # 초기화
            spec_res = {"passed": False, "off_target_count": 0, "details": "Skipped"}
            final_pass = False

            if thermo_res["passed"]:
                # -----------------------------------------------------
                # Step 2: Specificity QC (BLAST)
                # -----------------------------------------------------
                spec_res = self.spec_checker.check_amplicon(amp)
                
                # BLAST 통과 여부가 곧 최종 Specificity 결과가 됨
                final_pass = spec_res["passed"]
                
            else:
                spec_res["details"] = "Skipped (Thermo Failed)"
                final_pass = False

            # -----------------------------------------------------
            # 상태 업데이트
            # -----------------------------------------------------
            fail_reasons = []
            
            # 1. Thermo 실패 사유
            if not thermo_res["passed"]: 
                fail_reasons.append(thermo_res["fail_reason"])
            
            # 2. Specificity 실패 사유
            if thermo_res["passed"] and not spec_res["passed"]:
                # BlastChecker가 넣어둔 상세 정보 활용
                spec_details = amp.qc_details.get("specificity", {})
                
                # Signal Candidates(형광 발생) 개수로 사유 구체화
                # 1개면 정상(Intended), 0개면 Target 못 찾음, 2개 이상이면 Off-target
                cnt = spec_details.get("signal_candidates_count", 0)
                
                if cnt == 0:
                    fail_reasons.append("Spec(Target_Not_Found)")
                elif cnt > 1:
                    fail_reasons.append(f"Spec(Signal_OT:{cnt-1})")
                else:
                    fail_reasons.append("Specificity_Fail")
            

            # 3. QC Status 객체 생성
            amp.qc_status = AmpliconQCStatus(
                thermo=thermo_res,
                specificity=spec_res,
                ispcr=None, # 외부 isPCR 결과 없음
                overall_passed=final_pass
            )
            amp.is_qc_pass = final_pass
            
            # 전체 실패 사유 저장
            amp.qc_details["fail_reason"] = ", ".join(fail_reasons) if fail_reasons else ""
            
            processed_amplicons.append(amp)

        return processed_amplicons

    def summarize(self, amplicons: List[Amplicon]) -> pd.DataFrame:
        """
        엑셀 저장용 요약 데이터 생성
        """
        data = []
        for amp in amplicons:
            # 1. 기본 정보
            row = {
                "ID": amp.id,
                "Chrom": amp.reference_id,
                "Product_Size": amp.product_size,
                "QC_Result": "Pass" if amp.is_qc_pass else "Fail",
                "Fail_Reason": amp.qc_details.get("fail_reason", ""),
            }

            # 2. Specificity 상세 정보
            spec_info = amp.qc_details.get("specificity", {})
            if isinstance(spec_info, dict):
                # Count 정보가 가장 중요함
                row["Spec_Intended"] = "O" if spec_info.get("intended_found") else "X"
                row["Spec_Signal_OT"] = spec_info.get("signal_candidates_count", 0) # 전체 신호 발생 수
                row["Spec_Amp_OT"] = spec_info.get("amplification_only_count", 0)
                row["Spec_Total_Hits"] = spec_info.get("blast_hits_total", 0)
            else:
                row["Spec_Intended"] = "-"
                row["Spec_Signal_OT"] = 0
                row["Spec_Amp_OT"] = 0
                row["Spec_Total_Hits"] = 0

            # 3. Primer/Probe 서열, Tm, GC
            row.update({
                "Fwd_Seq": amp.forward.sequence,
                "Rev_Seq": amp.reverse.sequence,
                "Rev_Seq_RC": str(Seq(amp.reverse.sequence).reverse_complement()),
                "Fwd_Tm": round(amp.forward.tm, 2),
                "Rev_Tm": round(amp.reverse.tm, 2),
                "Fwd_GC": round(amp.forward.gc_percent, 2),
                "Rev_GC": round(amp.reverse.gc_percent, 2),
                # 서열 정보가 있다면 추가 (BLAST 단계에서 주입됨)
                "Reference_Seq": getattr(amp, "reference_sequence", ""),
                "Template_Seq": getattr(amp, "template_sequence", ""),
            })

            if amp.probe:
                p_seq = amp.probe.sequence
                p_gc = (p_seq.count('G') + p_seq.count('C')) / len(p_seq) * 100 if len(p_seq) > 0 else 0
                row.update({
                    "Probe_Seq": p_seq,
                    "Probe_Tm": round(amp.probe.tm, 2),
                    "Probe_GC": round(p_gc, 1)
                })
            else:
                row.update({"Probe_Seq": "N/A", "Probe_Tm": 0.0, "Probe_GC": 0.0})

            # 4. Thermo 상세 정보
            if amp.qc_status and amp.qc_status.thermo:
                t_data = amp.qc_status.thermo.get("data", {})
                row["HP_Fwd"] = round(t_data.get("fwd_hairpin_dg", 0.0), 2)
                row["HP_Rev"] = round(t_data.get("rev_hairpin_dg", 0.0), 2)
                row["HP_Probe"] = round(t_data.get("probe_hairpin_dg", 0.0), 2)
                row["Hetero_FR"] = round(t_data.get("hetero_fr_dg", 0.0), 2)
                row["Hetero_FP"] = round(t_data.get("hetero_fp_dg", 0.0), 2)
                row["Hetero_RP"] = round(t_data.get("hetero_rp_dg", 0.0), 2)

            data.append(row)
            
        df = pd.DataFrame(data)
        
        # 컬럼 순서 정리
        desired_order = [
            "ID", "Chrom", "QC_Result", "Fail_Reason", 
            "Spec_Intended", "Spec_Signal_OT", "Spec_Amp_OT", "Spec_Total_Hits",
            "Product_Size", 
            "Fwd_Tm", "Rev_Tm", "Probe_Tm", 
            "Fwd_GC", "Rev_GC", "Probe_GC",
            "HP_Fwd", "HP_Rev", "HP_Probe",
            "Hetero_FR", "Hetero_FP", "Hetero_RP",
            "Fwd_Seq", "Rev_Seq", "Rev_Seq_RC", "Probe_Seq", "Reference_Seq", "Template_Seq"
        ]
        
        final_cols = [c for c in desired_order if c in df.columns]
        remaining_cols = [c for c in df.columns if c not in final_cols]
        
        return df[final_cols + remaining_cols]