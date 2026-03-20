from pydantic import BaseModel, Field, model_validator
from typing import Dict, Any, List, Optional, Literal
from pcr.components.amplicon import Amplicon

# =================================================================
# [Input Schema]
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
        """
        [Adapter Method]
        Summary 객체와 Results 객체로 완벽히 분리된 계층형 JSON을 반환합니다.
        QC 항목이 유연하게 확장될 수 있도록 Workflow 패턴을 적용했습니다.
        """
        # 🔥 개선 2: 명확한 Summary 객체 생성
        summary = DesignSummary(
            status=self.status,
            total_count=self.total_count,
            passed_count=self.passed_count,
            failed_count=self.total_count - self.passed_count,
            error_msg=self.error_msg,
            log_messages=self.log_messages
        ).model_dump()

        if self.status != "success":
            return {"summary": summary, "results": []}

        results_list = []
        
        for rank, amp in enumerate(self.amplicons, start=1):
            
            # 🔥 개선 3: 유연한 QC 모듈 확장 (Workflow 패턴 지원)
            # 앞으로 어떤 QC 체커(ex: SNPChecker, DimerChecker)가 추가되든, 
            # 해당 체커가 amp.qc_metrics 딕셔너리에 결과를 담아두기만 하면 
            # 여기서 스키마 수정 없이 알아서 블랙박스 형태로 빨아들여 포장합니다.
            dynamic_qc_details = getattr(amp, "qc_metrics", {})
            
            # (하위 호환성) 기존에 존재하던 blast_stats를 동적 모듈 시스템으로 병합
            if hasattr(amp, "blast_stats") and getattr(amp, "blast_stats"):
                dynamic_qc_details["blast"] = getattr(amp, "blast_stats")

            # Amplicon 하나당 독립적인 객체 1개를 생성 (데이터의 원형 보존)
            result_item = {
                "rank": rank,
                "id": amp.id,
                
                # 핵심 판정 결과
                "qc_info": {
                    "is_pass": getattr(amp, "is_qc_pass", True),
                    "fail_reason": getattr(amp, "qc_log", "") if not getattr(amp, "is_qc_pass", True) else "PASS"
                },
                
                # 올리고(Oligo) 정보 그룹화
                "oligos": {
                    "forward": {
                        "sequence": amp.forward.sequence,
                        "tm": round(amp.forward.tm, 2),
                        "gc": round(getattr(amp.forward, "gc_percent", 0.0), 2)
                    },
                    "reverse": {
                        "sequence": amp.reverse.sequence,
                        "tm": round(amp.reverse.tm, 2),
                        "gc": round(getattr(amp.reverse, "gc_percent", 0.0), 2)
                    },
                    "probe": {
                        "sequence": amp.probe.sequence if amp.probe else "-",
                        "tm": round(amp.probe.tm, 2) if amp.probe else 0.0,
                        "gc": round(getattr(amp.probe, "gc_percent", 0.0), 2) if amp.probe else 0.0
                    }
                },
                
                # 앰플리콘 정보 그룹화
                "amplicon_info": {
                    "size": amp.product_size,
                    "tm": round(getattr(amp, "tm", 0.0), 2),
                    "gc": round(getattr(amp, "gc_percent", 0.0), 2),
                    "genomic_pos": getattr(amp, "genomic_pos", "Unknown"),
                    "alignment_text_block": "\n".join(getattr(amp, "alignment_visual", []))
                },
                
                # 🔥 무한히 확장 가능한 QC 상세 결과 컨테이너
                "qc_details": dynamic_qc_details
            }
            
            results_list.append(result_item)

        # 최종 반환 구조: Summary와 Results의 완벽한 2단 분리
        return {
            "summary": summary,
            "results": results_list
        }