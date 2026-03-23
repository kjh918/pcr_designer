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
# [Input Schema] 코어 엔진용 입력
# =================================================================
class BaseDesignInput(BaseModel):
    name: str
    template_sequence: str
    target_start: Optional[int] = None
    target_end: Optional[int] = None
    reference_sequence: Optional[str] = None 
    template_genomic_start: Optional[int] = None 
    template_genomic_end: Optional[int] = None 
    reference_name: str = "hg38"
    overrides: Dict[str, Any] = Field(default_factory=dict)

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
                if data.get('target_start') is None or data.get('target_end') is None:
                    raise ValueError("서열에 대괄호('[', ']')로 타겟을 지정하거나 좌표를 명시해야 합니다.")
                if data['target_start'] > data['target_end']:
                    raise ValueError(f"target_start({data['target_start']})는 target_end({data['target_end']})보다 클 수 없습니다.")
        return data

# =================================================================
# [Output Schema sub-models]
# 프론트엔드 응답을 위한 명확한 객체 분리 (Summary / Results)
# =================================================================
class DesignSummary(BaseModel):
    status: Literal["success", "fail", "error"]
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
    status: Literal["success", "fail", "error"] 
    log_messages: List[str] = []
    error_msg: Optional[str] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)

    @property
    def passed_amplicons(self) -> List[Amplicon]:
        return [amp for amp in self.amplicons if getattr(amp, "is_qc_pass", False)]

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
            return {
                "status": self.status,
                "summary": summary,
                "metadata": self.metadata,
                "results": []
            }

        results_list = []
        for rank, amp in enumerate(self.amplicons, start=1):
            dynamic_qc_details = getattr(amp, "qc_metrics", {})
            if hasattr(amp, "blast_stats") and getattr(amp, "blast_stats"):
                dynamic_qc_details["blast"] = getattr(amp, "blast_stats")

            # 🔥 [수정] Heterodimer는 두 객체 사이의 관계값이므로 여전히 Metrics나 QC Status에서 가져옵니다.
            t_data = {}
            qc_status = getattr(amp, "qc_status", None)
            if qc_status and hasattr(qc_status, "modules") and "thermo" in qc_status.modules:
                t_data = qc_status.modules["thermo"].metrics

            # 헬퍼 함수: 객체의 속성을 안전하게 읽어옴 (getattr 사용)
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
                "qc_info": {
                    "is_pass": getattr(amp, "is_qc_pass", True),
                    "fail_reason": getattr(amp, "qc_log", "") if not getattr(amp, "is_qc_pass", True) else "PASS"
                },
                "oligos": {
                    "forward": get_oligo_meta(amp.forward),
                    "reverse": get_oligo_meta(amp.reverse),
                    "probe": get_oligo_meta(amp.probe) if amp.probe else {
                        "sequence": "-", "tm": 0, "gc": 0, "cpg_count": 0, "hairpin_dg": 0, "homodimer_dg": 0
                    },
                    "heterodimer": {
                        "fr_dg": round(t_data.get("hetero_fr_dg", 0.0), 2),
                        "fp_dg": round(t_data.get("hetero_fp_dg", 0.0), 2),
                        "rp_dg": round(t_data.get("hetero_rp_dg", 0.0), 2)
                    }
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

        return {
            "status": self.status,
            "metadata": self.metadata,
            "summary": summary,
            "results": results_list
        }