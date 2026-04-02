from typing import Optional, Dict, Any
from pydantic import BaseModel, Field, model_validator
from pcr.designers.base.schema import BaseDesignInput, BaseDesignOutput, DesignSummary

class ASPCRDesignInput(BaseDesignInput):
    design_name: str = "AS-PCR_Design"
    sequence: str = "" 
    reference_genome: str = "hg38"
    top_k: int = 5
    
    template_genomic_start: int = 0
    template_genomic_end: int = 0

    # 🔥 웹에서 넘어오는 AS-PCR 전용 설정
    fixed_prime: str = "forward"
    mismatch_pos: int = 3
    mismatch_intensity: str = "strong"

    amplicon: Dict[str, Any] = Field(default_factory=dict)
    primer: Dict[str, Any] = Field(default_factory=dict)
    qc_criteria: Dict[str, Any] = Field(default_factory=dict)

    templates: Dict[str, str] = Field(default_factory=dict)

    @model_validator(mode='before')
    @classmethod
    def parse_brackets_and_validate(cls, data: Any) -> Any:
        if isinstance(data, dict):
            data['name'] = data.get('design_name', 'AS-PCR_Design')
            data['reference_name'] = data.get('reference_genome', 'hg38')
            if 'sequence' in data:
                data['template_sequence'] = data['sequence']
            data['target_start'] = data.get('target_start', 0)
            data['target_end'] = data.get('target_end', 0)
        return data

    def to_core_pcr_params(self) -> Dict[str, Any]:
        """🔥 프론트엔드 파라미터를 코어 엔진용 파라미터로 완벽 맵핑"""
        return {
            "primer_kwargs": {
                "n_candidates": 5,
                "min_amplicon_length": self.amplicon.get("min_length", 60),
                "max_amplicon_length": self.amplicon.get("max_length", 150),
                "min_length": self.primer.get("min_length", 15),
                "opt_length": self.primer.get("opt_length", 22),
                "max_length": self.primer.get("max_length", 30),
                "min_tm": self.primer.get("min_tm", 52.0),
                "opt_tm": self.primer.get("opt_tm", 58.0),
                "max_tm": self.primer.get("max_tm", 65.0),
                "min_gc": self.primer.get("min_gc", 35.0),
                "opt_gc": self.primer.get("opt_gc", 50.0),
                "max_gc": self.primer.get("max_gc", 65.0),
            }
        }

    def to_core_qc_overrides(self) -> Dict[str, Any]:
        """🔥 예시로 주신 코드처럼 모든 QC 옵션과 Alias를 빠짐없이 매핑합니다."""
        thermo = self.qc_criteria.get("thermodynamics", {})
        blast = self.qc_criteria.get("blast", {})
        amp = self.qc_criteria.get("amplicon", {})
        oligo = self.qc_criteria.get("oligo", {})
        
        overrides = {
            # Thermodynamics
            "hairpin_min_dg": thermo.get("hairpin_min_dg", -5.0),
            "homodimer_min_dg": thermo.get("homodimer_min_dg", -6.0),
            "heterodimer_min_dg": thermo.get("heterodimer_min_dg", -6.0),
            
            # BLAST (Alias 포함)
            "min_identity": blast.get("min_identity", 90.0),
            "blast_identity_threshold": blast.get("min_identity", 90.0),
            "min_hit_length": blast.get("min_hit_length", 13),
            "blast_max_alignments": blast.get("max_alignments", 50),
            "max_alignments": blast.get("max_alignments", 50),
            
            # Amplicon Size (Alias 포함)
            "min_amp_size": amp.get("min_size", 50),
            "min_amp_len": amp.get("min_size", 50),
            "max_amp_size": amp.get("max_size", 300),
            "max_amp_len": amp.get("max_size", 300),
            
            "use_ispcr_check": amp.get("use_ispcr", False),
            
            # Primer Constraints
            "primer": {
                "max_diff_tm": oligo.get("max_tm_diff", 5.0)
            },
            
            # 🔥 AS-PCR 스크립트 단에서 접근할 수 있도록 포장
            "as_pcr": {
                "fixed_prime": self.fixed_prime,
                "mismatch_pos": self.mismatch_pos,
                "mismatch_intensity": self.mismatch_intensity
            }
        }
        return overrides

class ASPCRDesignOutput(BaseDesignOutput):
    def to_frontend_dict(self) -> Dict[str, Any]:
        summary = DesignSummary(
            status=self.status, total_count=self.total_count,
            passed_count=self.passed_count, failed_count=self.total_count - self.passed_count,
            error_msg=self.error_msg, log_messages=self.log_messages
        ).model_dump()

        if self.status != "success":
            return {"status": self.status, "summary": summary, "metadata": self.metadata, "results": []}

        sets_dict = {}
        for amp in self.amplicons:
            set_id = getattr(amp, "set_id", "UnknownSet")
            if set_id not in sets_dict: sets_dict[set_id] = []
            sets_dict[set_id].append(amp)

        results_list = []
        for rank, (set_id, amps_in_set) in enumerate(sets_dict.items(), start=1):
            set_qc_pass = all(getattr(amp, "is_qc_pass", False) for amp in amps_in_set)
            set_data = {"rank": rank, "set_id": set_id, "set_qc_pass": set_qc_pass, "alleles": {}}
            
            for amp in amps_in_set:
                allele_type = getattr(amp, "allele_type", "unknown")
                dynamic_qc_details = getattr(amp, "qc_metrics", {}).copy()
                blast_stats = getattr(amp, "blast_stats", {})
                if blast_stats: dynamic_qc_details["blast"] = blast_stats

                t_data = {}
                qc_status = getattr(amp, "qc_status", None)
                if qc_status and hasattr(qc_status, "modules") and "thermo" in qc_status.modules:
                    t_data = qc_status.modules["thermo"].metrics

                if qc_status is not None:
                    final_is_pass = qc_status.is_pass
                    final_fail_reason = " | ".join(qc_status.fail_reasons) if not final_is_pass else "PASS"
                else:
                    final_is_pass = getattr(amp, "is_qc_pass", True)
                    final_fail_reason = getattr(amp, "qc_log", "PASS") if not final_is_pass else "PASS"

                if getattr(amp, "is_qc_pass", None) is False:
                    final_is_pass = False
                    final_fail_reason = getattr(amp, "qc_log", final_fail_reason)

                def get_oligo_meta(obj):
                    if not obj: return {"tm": 0.0, "gc": 0.0, "cpg": 0, "hp": 0.0, "hd": 0.0}
                    return {
                        "sequence": getattr(obj, "sequence", "-"),
                        "tm": round(getattr(obj, "tm", 0.0), 2),
                        "gc": round(getattr(obj, "gc_percent", 0.0), 2),
                        "cpg_count": int(getattr(obj, "cpg_count", 0)),
                        "hairpin_dg": round(getattr(obj, "hairpin_dg", 0.0), 2),
                        "homodimer_dg": round(getattr(obj, "homodimer_dg", 0.0), 2)
                    }

                aln_text = "\n".join(getattr(amp, "alignment_visual", []))
                
                if blast_stats and "_OffTarget_" not in amp.id and "_Noise_" not in amp.id:
                    ot_sigs = blast_stats.get("off_target_signals", [])
                    noise_sigs = blast_stats.get("amplification_only", [])
                    
                    if ot_sigs or noise_sigs:
                        extra = ["\n", "=" * 65, "⚠️ NON-SPECIFIC BINDINGS (OFF-TARGET / NOISE)", "=" * 65]
                        for ot in ot_sigs:
                            extra.extend(["\n[Off-Target : Probe Binds Here]", ot.get("unified_text_block", "")])
                        for noise in noise_sigs:
                            extra.extend(["\n[Noise : Amplification Only]", noise.get("unified_text_block", "")])
                        
                        aln_text += "\n".join(extra)

                set_data["alleles"][allele_type] = {
                    "id": amp.id,
                    "metrics": {"pair_penalty": round(getattr(amp, "pair_penalty", 0), 3)},
                    "qc_info": {"is_pass": final_is_pass, "fail_reason": final_fail_reason},
                    "oligos": {
                        "forward": get_oligo_meta(amp.forward),
                        "reverse": get_oligo_meta(amp.reverse),
                        "probe": get_oligo_meta(amp.probe) if amp.probe else {"sequence": "-", "tm": 0, "gc": 0, "cpg_count": 0, "hairpin_dg": 0, "homodimer_dg": 0},
                        "heterodimer": {
                            "fr_dg": round(t_data.get("hetero_fr_dg", 0.0), 2), 
                            "fp_dg": round(t_data.get("hetero_fp_dg", 0.0), 2), 
                            "rp_dg": round(t_data.get("hetero_rp_dg", 0.0), 2)
                        }
                    },
                    "amplicon_info": {
                        "sequence": getattr(amp, "sequence", "-"),
                        "size": amp.product_size,
                        "tm": round(getattr(amp, "tm", 0.0), 2),
                        "gc": round(getattr(amp, "gc_percent", 0.0), 2),
                        "genomic_pos": getattr(amp, "genomic_pos", "Unknown"),
                        "alignment_text_block": aln_text
                    },
                    "qc_details": dynamic_qc_details
                }
            
            results_list.append(set_data)

        summary["total_count"] = len(results_list)
        summary["passed_count"] = sum(1 for s in results_list if s["set_qc_pass"])
        summary["failed_count"] = summary["total_count"] - summary["passed_count"]

        return {"status": self.status, "metadata": self.metadata, "summary": summary, "results": results_list}