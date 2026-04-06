from typing import Optional, Dict, Any, List, Union
from pydantic import BaseModel, Field, model_validator
from pcr.designers.base.schema import BaseDesignInput, BaseDesignOutput, DesignSummary

# =================================================================
# [Layer 1] 하위 파라미터 모델 정의 (AttributeError 방지용)
# =================================================================

class MspcrPrimerParams(BaseModel):
    min_length: int = 15
    opt_length: int = 18
    max_length: int = 25
    min_tm: float = 45.0
    opt_tm: float = 55.0
    max_tm: float = 65.0
    min_gc: float = 25.0
    opt_gc: float = 45.0
    max_gc: float = 85.0

class MspcrAmpliconParams(BaseModel):
    min_length: int = 80
    max_length: int = 200

class MspcrQcThermodynamics(BaseModel):
    hairpin_min_dg: float = -6.0
    homodimer_min_dg: float = -6.0
    heterodimer_min_dg: float = -6.0

class MspcrQcBlast(BaseModel):
    min_identity: float = 80.0
    min_hit_length: int = 10
    max_alignments: int = 10

class MspcrQcAmplicon(BaseModel):
    min_size: int = 50
    max_size: int = 300
    use_ispcr: bool = False

class MspcrQcOligo(BaseModel):
    max_tm_diff: float = 3.0
    window_size_3prime: int = 3
    min_cpg_count: int = 1

class MspcrQcCriteriaModel(BaseModel):
    """라우터에서 호출하는 to_core_overrides 메서드를 보유한 모델"""
    thermodynamics: MspcrQcThermodynamics = Field(default_factory=MspcrQcThermodynamics)
    blast: MspcrQcBlast = Field(default_factory=MspcrQcBlast)
    amplicon: MspcrQcAmplicon = Field(default_factory=MspcrQcAmplicon)
    oligo: MspcrQcOligo = Field(default_factory=MspcrQcOligo)

    def to_core_overrides(self) -> Dict[str, Any]:
        """
        계층형 모델의 내부 속성 경로를 지정하여 딕셔너리로 변환합니다.
        🔥 BLAST 섹션은 구조상 포함하되, 실행은 False로 고정합니다.
        """
        return {
            # Thermodynamics 섹션
            "hairpin_min_dg": self.thermodynamics.hairpin_min_dg,
            "homodimer_min_dg": self.thermodynamics.homodimer_min_dg,
            "heterodimer_min_dg": self.thermodynamics.heterodimer_min_dg,
            
            # BLAST 섹션 (구조 유지를 위해 값은 넣되, run_blast 플래그를 False로 고정)
            "min_identity": self.blast.min_identity,
            "blast_identity_threshold": self.blast.min_identity, 
            "min_hit_length": self.blast.min_hit_length,
            "blast_max_alignments": self.blast.max_alignments,
            "max_alignments": self.blast.max_alignments,
            "run_blast": False,  # 🔥 MS-PCR 전용 고정 설정
            
            # Amplicon 섹션
            "min_amp_size": self.amplicon.min_size,
            "max_amp_size": self.amplicon.max_size,
            "use_ispcr_check": self.amplicon.use_ispcr,
            
            # Oligo 섹션
            "primer": {"max_diff_tm": self.oligo.max_tm_diff}
        }

# =================================================================
# [Layer 2] 메인 입력 모델 (MS-PCR 전용)
# =================================================================

class MSPCRDesignInput(BaseDesignInput):
    """
    MS-PCR 전용 입력 객체.
    BaseDesignInput을 상속받으며, 계층형 파라미터를 모델로 수신합니다.
    """
    design_name: str = "MS-PCR_Design"
    sequence: str = "" 
    reference_genome: str = "hg38"
    top_k: int = 5
    
    # MS-PCR 전용 설정
    window_size_3prime: int = Field(default=3)
    min_cpg_count: int = Field(default=1)

    # 단순 Dict가 아닌 위에서 정의한 Model 클래스 사용
    amplicon: MspcrAmpliconParams = Field(default_factory=MspcrAmpliconParams)
    primer: MspcrPrimerParams = Field(default_factory=MspcrPrimerParams)
    qc_criteria: MspcrQcCriteriaModel = Field(default_factory=MspcrQcCriteriaModel)

    templates: Dict[str, str] = Field(default_factory=dict)
    target_cpg_indices: List[int] = Field(default_factory=list)

    @model_validator(mode='before')
    @classmethod
    def parse_mspcr_inputs(cls, data: Any) -> Any:
        if isinstance(data, dict):
            # API 필드 매핑
            data['name'] = data.get('design_name', 'MS-PCR_Design')
            data['reference_name'] = data.get('reference_genome', 'hg38')
            if 'sequence' in data:
                data['template_sequence'] = data['sequence']
            # 부모 클래스 필수 좌표값 기본값 설정
            data.setdefault('target_start', 0)
            data.setdefault('target_end', 0)
        return data

    def to_core_pcr_params(self) -> Dict[str, Any]:
        """프론트엔드 파라미터를 Primer3 엔진 규격으로 매핑"""
        return {
            "primer_kwargs": {
                "n_candidates": 30,
                "min_amplicon_length": self.amplicon.min_length,
                "max_amplicon_length": self.amplicon.max_length,
                "min_length": self.primer.min_length,
                "opt_length": self.primer.opt_length,
                "max_length": self.primer.max_length,
                "min_tm": self.primer.min_tm,
                "opt_tm": self.primer.opt_tm,
                "max_tm": self.primer.max_tm,
                "min_gc": self.primer.min_gc,
                "opt_gc": self.primer.opt_gc,
                "max_gc": self.primer.max_gc,
            }
        }

    def to_core_qc_overrides(self) -> Dict[str, Any]:
        """
        일반 QC 옵션에 MS-PCR 전용 파라미터(window_size_3prime 등)를 
        'ms_pcr' 키에 담아 포함합니다.
        """
        overrides = self.qc_criteria.to_core_overrides()
        
        # MS-PCR 전용 알고리즘을 위한 설정값 주입
        overrides["ms_pcr"] = {
            "window_size_3prime": self.window_size_3prime,
            "min_cpg_count": self.min_cpg_count
        }
        return overrides

# =================================================================
# [Layer 3] 결과 출력 모델
# =================================================================

class MspcrDesignOutput(BaseDesignOutput):
    """MS-PCR 결과를 M/U 세트 구조로 변환하여 프론트엔드에 전달"""
    def to_frontend_dict(self) -> Dict[str, Any]:
        summary = DesignSummary(
            status=self.status, total_count=self.total_count,
            passed_count=self.passed_count, failed_count=self.total_count - self.passed_count,
            error_msg=self.error_msg, log_messages=self.log_messages
        ).model_dump()

        if self.status != "success":
            return {"status": self.status, "summary": summary, "metadata": self.metadata, "results": []}

        # 앰플리콘들을 Set ID별로 그룹화 (M/U 세트 구성)
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
                
                if blast_stats and "_OffTarget_" not in amp.id:
                    ot_sigs = blast_stats.get("off_target_signals", [])
                    noise_sigs = blast_stats.get("amplification_only", [])
                    if ot_sigs or noise_sigs:
                        extra = ["\n", "=" * 65, "⚠️ NON-SPECIFIC BINDINGS (OFF-TARGET / NOISE)", "=" * 65]
                        for ot in ot_sigs:
                            extra.extend(["\n[Off-Target : Potential Binding]", ot.get("unified_text_block", "")])
                        aln_text += "\n".join(extra)

                set_data["alleles"][allele_type] = {
                    "id": amp.id,
                    "metrics": {"pair_penalty": round(getattr(amp, "pair_penalty", 0), 3)},
                    "qc_info": {"is_pass": final_is_pass, "fail_reason": final_fail_reason},
                    "oligos": {
                        "forward": get_oligo_meta(amp.forward),
                        "reverse": get_oligo_meta(amp.reverse),
                        "heterodimer": {
                            "fr_dg": round(t_data.get("hetero_fr_dg", 0.0), 2)
                        }
                    },
                    "amplicon_info": {
                        "sequence": getattr(amp, "sequence", "-"),
                        "size": amp.product_size,
                        "tm": round(getattr(amp, "tm", 0.0), 2),
                        "gc": round(getattr(amp, "gc_percent", 0.0), 2),
                        "cpg_count": int(getattr(amp, "cpg_count", 0)),
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