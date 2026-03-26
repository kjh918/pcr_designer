import copy
from pydantic import BaseModel, Field, model_validator
from typing import Dict, Any, List, Optional, Literal
from pcr.components.amplicon import Amplicon

# =================================================================
# [API Request Schema] 프론트엔드 데이터 수신용 뼈대
# =================================================================
class SequencesModel(BaseModel):
    forward: str
    reverse: str
    probe: str = ""
    template: str = ""

class QCEvalInput(BaseModel):
    project_name: str = "QC_Project"
    sequences: SequencesModel
    reference_genome: str = "none"
    qc_criteria: Dict[str, Any] = Field(default_factory=dict)


# =================================================================
# [Input Schema] 코어 엔진용 입력 (Base)
# =================================================================
class BaseDesignInput(BaseModel):
    """
    모든 기법(qPCR, MS-PCR, AS-PCR 등)이 공통으로 사용하는 순수 뼈대입니다.
    기법별 특수 파라미터는 여기 두지 않고, 상속받는 자식 클래스에 정의합니다.
    """
    name: str
    template_sequence: str
    target_start: Optional[int] = None
    target_end: Optional[int] = None
    reference_sequence: Optional[str] = None 
    template_genomic_start: Optional[int] = None 
    template_genomic_end: Optional[int] = None 
    reference_name: str = "hg38"
    overrides: Dict[str, Any] = Field(default_factory=dict)

    # Factory 실행에 필수적인 Config
    config: Optional[Any] = Field(default=None, description="PCRFactory에서 주입되는 PipelineConfig 객체 (필수)")
    
    # Factory 공통 파싱 결과물 (대괄호 파싱)
    target_indices: List[int] = Field(default_factory=list, description="대괄호 파싱으로 추출된 타겟(SNP/CpG) 인덱스 목록")

    @model_validator(mode='before')
    @classmethod
    def parse_brackets_and_validate(cls, data: Any) -> Any:
        if isinstance(data, dict):
            seq = data.get('template_sequence', '')
            if '[' in seq and ']' in seq:
                start_idx = seq.find('[')
                end_idx = seq.find(']') - 1
                clean_seq = seq.replace('[', '').replace(']', '')
                if len(clean_seq) > 2000:
                    raise ValueError(f"Template Sequence는 2K 이하입니다.")
                data['template_sequence'] = clean_seq
                data['target_start'] = start_idx
                data['target_end'] = end_idx
            else:
                if data.get('target_start') is not None and data.get('target_end') is not None:
                    if data['target_start'] > data['target_end']:
                        raise ValueError(f"target_start({data['target_start']})는 target_end({data['target_end']})보다 클 수 없습니다.")
        return data

# =================================================================
# [Output Schema sub-models]
# =================================================================
class DesignSummary(BaseModel):
    status: str  # 🔥 Literal 강제 검증을 해제하고 str로 유연하게 변경
    total_count: int = 0
    passed_count: int = 0
    failed_count: int = 0
    error_msg: Optional[str] = None
    log_messages: List[str] = []

# =================================================================
# [Main Output Schema]
# =================================================================
class BaseDesignOutput(BaseModel):
    amplicons: List[Amplicon] = []
    status: str  # 🔥 Literal 강제 검증 해제 ("no_probes_found" 등 자유로운 상태값 허용)
    log_messages: List[str] = []
    error_msg: Optional[str] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)

    @property
    def passed_amplicons(self) -> List[Amplicon]:
        passed = []
        for amp in self.amplicons:
            qs = getattr(amp, "qc_status", None)
            if qs is not None:
                if qs.is_pass: 
                    passed.append(amp)
            else:
                if getattr(amp, "is_qc_pass", False): 
                    passed.append(amp)
        return passed

    @property
    def total_count(self) -> int:
        return len(self.amplicons)

    @property
    def passed_count(self) -> int:
        return len(self.passed_amplicons)

    def to_frontend_dict(self) -> Dict[str, Any]:
        summary = DesignSummary(
            status=self.status,
            total_count=self.total_count,
            passed_count=self.passed_count,
            failed_count=self.total_count - self.passed_count,
            error_msg=self.error_msg,
            log_messages=self.log_messages
        ).model_dump()

        if self.status != "success":
            return {"status": self.status, "summary": summary, "metadata": self.metadata, "results": []}

        results_list = []
        for rank, amp in enumerate(self.amplicons, start=1):
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

            result_item = {
                "rank": rank,
                "id": amp.id,
                "qc_info": {"is_pass": final_is_pass, "fail_reason": final_fail_reason},
                "oligos": {
                    "forward": get_oligo_meta(amp.forward),
                    "reverse": get_oligo_meta(amp.reverse),
                    "probe": get_oligo_meta(amp.probe) if amp.probe else {"sequence": "-", "tm": 0, "gc": 0, "cpg_count": 0, "hairpin_dg": 0, "homodimer_dg": 0},
                    "heterodimer": {"fr_dg": round(t_data.get("hetero_fr_dg", 0.0), 2), "fp_dg": round(t_data.get("hetero_fp_dg", 0.0), 2), "rp_dg": round(t_data.get("hetero_rp_dg", 0.0), 2)}
                },
                "amplicon_info": {
                    "size": amp.product_size,
                    "tm": round(getattr(amp, "tm", 0.0), 2),
                    "gc": round(getattr(amp, "gc_percent", 0.0), 2),
                    "genomic_pos": getattr(amp, "genomic_pos", "Unknown"),
                    "alignment_text_block": "\n".join(getattr(amp, "alignment_visual", []))
                },
                "qc_details": dynamic_qc_details
            }
            results_list.append(result_item)

        return {"status": self.status, "metadata": self.metadata, "summary": summary, "results": results_list}