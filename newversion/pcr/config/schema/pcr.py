from typing import Dict, Any, Optional, Literal
from pydantic import BaseModel, Field

# =========================================================================
# 1. 공통 설정을 위한 기본 클래스 (Primer/Probe 공통)
# =========================================================================
class BaseDesignKwargs(BaseModel):
    """
    Primer와 Probe가 공통으로 가지는 속성 (길이, Tm, GC, 개수)
    기본값을 지정하여 YAML에 누락되어도 에러가 나지 않도록 함.
    """
    # Length
    opt_length: int = Field(20, description="Optimal Length")
    min_length: int = Field(18, description="Min Length")
    max_length: int = Field(25, description="Max Length")
    
    # Tm
    opt_tm: float = Field(50.0, description="Optimal Tm")
    min_tm: float = Field(55.0, description="Min Tm")
    max_tm: float = Field(60.0, description="Max Tm")
    
    # GC
    opt_gc: float = Field(40.0, description="Optimal GC %")
    min_gc: float = Field(50.0, description="Min GC %")
    max_gc: float = Field(60.0, description="Max GC %")
    
    # Candidates (YAML의 PRIMER_NUM_RETURN 등과 매핑)
    n_candidates: int = Field(10, alias="PRIMER_NUM_RETURN")

    def to_primer3_args(self, prefix: str = "PRIMER_") -> Dict[str, Any]:
        """
        Pydantic 필드를 Primer3 옵션 키로 변환
        prefix: "PRIMER_" (Primer용) 또는 "PRIMER_INTERNAL_" (Probe용)
        """
        return {
            f"{prefix}OPT_SIZE": self.opt_length,
            f"{prefix}MIN_SIZE": self.min_length,
            f"{prefix}MAX_SIZE": self.max_length,
            f"{prefix}OPT_TM": self.opt_tm,
            f"{prefix}MIN_TM": self.min_tm,
            f"{prefix}MAX_TM": self.max_tm,
            f"{prefix}OPT_GC_PERCENT": self.opt_gc,
            f"{prefix}MIN_GC": self.min_gc,
            f"{prefix}MAX_GC": self.max_gc,
            # Probe인 경우 호출 측에서 Key를 PRIMER_INTERNAL_NUM_RETURN으로 덮어써야 할 수 있음
            f"{prefix}NUM_RETURN": self.n_candidates 
        }

# =========================================================================
# 2. Primer 설정
# =========================================================================
class PrimerKwargs(BaseDesignKwargs):
    """Primer 전용 설정 (Product Size 등)"""
    
    # Primer3 Product Size Range
    min_amplicon_length: int = Field(80, description="Min Product Size")
    max_amplicon_length: int = Field(150, description="Max Product Size")
    
    # 추가적인 Primer3 Global Args (YAML에서 자유롭게 추가 가능)
    primer3_global_args: Dict[str, Any] = Field(default_factory=dict)

    def to_global_args(self) -> Dict[str, Any]:
        """Primer 설계를 위한 전체 Global Args 생성 (Prefix: PRIMER_)"""
        # 1. Base 속성 변환
        args = self.to_primer3_args(prefix="PRIMER_")
        
        # 2. Primer 전용 속성 추가
        args.update({
            "PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
            "PRIMER_NUM_RETURN": self.n_candidates,
        })
        
        # 3. 기타 사용자 정의 Args 병합
        if self.primer3_global_args:
            args.update(self.primer3_global_args)
            
        return args

# =========================================================================
# 3. Probe 설정
# =========================================================================
class ProbeKwargs(BaseDesignKwargs):
    """Probe 전용 설정 (Tm이 더 높음, Python 필터링 옵션 포함)"""
    
    # [Override Defaults] Probe는 보통 Primer보다 Tm이 높고 길이가 긺
    opt_length: int = 25
    min_length: int = 20
    max_length: int = 30
    
    opt_tm: float = 65.0
    min_tm: float = 60.0
    max_tm: float = 70.0
    
    min_gc: float = 40.0
    max_gc: float = 60.0
        
    n_candidates: int = Field(100, alias="PRIMER_NUM_RETURN") # Base와 이름 통일
        
    # [Tm Diff Logic]
    min_primer_probe_tm_diff: float = 5.0
    max_primer_probe_tm_diff: float = 10.0

    # ✅ [NEW] Python Logic용 Custom Filter Options
    max_probe_poly_g: int = Field(
        4, description="Probe 내 연속된 G의 최대 허용 개수 (예: 4면 GGGG 허용)"
    )
    max_probe_3_end_gc: int = Field(
        3, description="Probe 3' 말단 5bp 내 G/C의 최대 개수"
    )

    def to_global_args(self) -> Dict[str, Any]:
        """Probe 설계를 위한 Global Args 생성 (Prefix: PRIMER_INTERNAL_)"""
        # 1. Base 매핑 (PRIMER_INTERNAL_ 접두어 사용)
        args = self.to_primer3_args(prefix="PRIMER_INTERNAL_")
        
        # 2. Key 보정 (Primer3는 Probe 개수를 PRIMER_INTERNAL_NUM_RETURN으로 받음)
        args["PRIMER_INTERNAL_NUM_RETURN"] = self.n_candidates
        
        # 참고: max_probe_poly_g 등은 Primer3 옵션이 아니므로 여기엔 포함 안 됨 (Python 로직에서 사용)
        return args

# =========================================================================
# 4. 기타 설정 (Bisulfite, AS-PCR)
# =========================================================================
class BisulfiteKwargs(BaseModel):
    run: bool = False
    cpg_default: str = "methyl" # unmethyl, both

class ASPCRKwargs(BaseModel):
    target_allele: Literal["ref", "alt"] = Field("alt", description="Target allele to amplify")
        
    # Artificial Mismatch (AM) 설정
    use_artificial_mismatch: bool = Field(
        False, description="3' 말단 근처에 인위적 불일치 도입 여부"
    )
    
    # 위치: -2 (penultimate), -3 (antepenultimate)
    am_position: Literal["-2", "-3"] = Field("-3")
        
    # 강도: Strong/Weak/Auto
    am_strength: Literal["strong", "weak", "auto"] = Field("auto")

    def get_position_index(self) -> int:
        """Python Slicing용 인덱스 반환 (-2 or -3)"""
        return int(self.am_position)

# =========================================================================
# 5. 통합 PCR Params (Root)
# =========================================================================
class PCRParams(BaseModel):
    # 필수
    primer_kwargs: PrimerKwargs
    
    # 선택 사항 (Optional + Default None 처리로 에러 방지)
    probe_kwargs: Optional[ProbeKwargs] = None
    bisulfite_kwargs: Optional[BisulfiteKwargs] = None
    as_pcr_kwargs: Optional[ASPCRKwargs] = None