from pydantic import BaseModel
from typing import Dict, Any, Optional

class PrimerKwargs(BaseModel):
    n_candidates: int = 100
    min_amplicon_length: int = 80
    max_amplicon_length: int = 200
    opt_length: int = 20
    min_length: int = 18
    max_length: int = 27
    opt_tm: float = 60.0
    min_tm: float = 55.0
    max_tm: float = 65.0
    opt_gc: float = 45.0
    min_gc: float = 25.0
    max_gc: float = 85.0

    def to_global_args(self, *args, **kwargs) -> Dict[str, Any]:
        return {
            "PRIMER_NUM_RETURN": self.n_candidates,
            "PRIMER_OPT_SIZE": self.opt_length,
            "PRIMER_MIN_SIZE": self.min_length,
            "PRIMER_MAX_SIZE": self.max_length,
            "PRIMER_OPT_TM": self.opt_tm,
            "PRIMER_MIN_TM": self.min_tm,
            "PRIMER_MAX_TM": self.max_tm,
            "PRIMER_OPT_GC_PERCENT": self.opt_gc,
            "PRIMER_MIN_GC": self.min_gc,
            "PRIMER_MAX_GC": self.max_gc,
            "PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]]
        }
    
class ProbeKwargs(BaseModel):
    n_candidates: int = 100
    n_primers_per_probe: int = 10 
    probe_gap: int = 3            
    opt_length: int = 25
    min_length: int = 20
    max_length: int = 30
    
    # 동적 Tm 계산용 Diff 변수
    opt_tm_diff: float = 7.0
    min_tm_diff: float = 5.0
    max_tm_diff: float = 10.0
    
    opt_gc: float = 50.0
    min_gc: float = 40.0
    max_gc: float = 60.0
    
    # 💡 [핵심 수정] 어떤 값이 들어오든 유연하게 Primer의 opt_tm을 뽑아냅니다.
    def to_global_args(self, primer_arg: Any = None, **kwargs) -> Dict[str, Any]:
        # 1. 기본 프라이머 Tm 기준값
        base_tm = 60.0 

        # 2. 키워드로 명시해서 보낸 경우 (예: primer_opt_tm=60.0)
        if "primer_opt_tm" in kwargs:
            base_tm = float(kwargs["primer_opt_tm"])
            
        # 3. BaseDesigner에서 PrimerKwargs 객체를 통째로 던진 경우
        elif primer_arg is not None:
            if hasattr(primer_arg, "opt_tm"):
                base_tm = float(primer_arg.opt_tm) # 객체 안에서 opt_tm만 쏙 빼옵니다.
            elif isinstance(primer_arg, (int, float)):
                base_tm = float(primer_arg)
                
        return {
            "PRIMER_INTERNAL_OPT_SIZE": self.opt_length,
            "PRIMER_INTERNAL_MIN_SIZE": self.min_length,
            "PRIMER_INTERNAL_MAX_SIZE": self.max_length,
            
            # 🔥 객체 에러 없이 완벽하게 실수(Float) 끼리 덧셈이 이루어짐
            "PRIMER_INTERNAL_OPT_TM": base_tm + self.opt_tm_diff,
            "PRIMER_INTERNAL_MIN_TM": base_tm + self.min_tm_diff,
            "PRIMER_INTERNAL_MAX_TM": base_tm + self.max_tm_diff,
            
            "PRIMER_INTERNAL_OPT_GC_PERCENT": self.opt_gc,
            "PRIMER_INTERNAL_MIN_GC": self.min_gc,
            "PRIMER_INTERNAL_MAX_GC": self.max_gc,
        }

class PCRParams(BaseModel):
    primer_kwargs: PrimerKwargs
    probe_kwargs: Optional[ProbeKwargs] = None