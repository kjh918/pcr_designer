from pydantic import BaseModel, Field
from typing import Dict, Any, Optional, Literal

# -------------------------------------------------------------------------
# 1. 공통 설정을 위한 기본 클래스
# -------------------------------------------------------------------------
class BaseDesignKwargs(BaseModel):
    opt_length: int = 20
    min_length: int = 18
    max_length: int = 27
    
    opt_tm: float = 60.0
    min_tm: float = 55.0
    max_tm: float = 65.0
    
    opt_gc: float = 45.0
    min_gc: float = 25.0
    max_gc: float = 85.0
    
    n_primers: int = 100

    def to_primer3_args(self, prefix: str = "PRIMER_") -> Dict[str, Any]:
        """
        YAML의 friendly name을 Primer3 옵션명으로 변환
        prefix를 인자로 받아 Primer("PRIMER_")와 Probe("PRIMER_INTERNAL_") 모두 대응
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
            f"{prefix}NUM_RETURN": self.n_primers,
        }

# -------------------------------------------------------------------------
# 2. Primer 설정
# -------------------------------------------------------------------------
class PrimerKwargs(BaseDesignKwargs):
    min_amplicon_length: int = 80
    max_amplicon_length: int = 120
    n_primers: int = 100
    primer3_global_args: Dict[str, Any] = Field(default_factory=dict)

    def to_global_args(self) -> Dict[str, Any]:
        """Primer 설계를 위한 전체 Global Args 생성 (Prefix: PRIMER_)"""
        args = self.to_primer3_args(prefix="PRIMER_")
        args.update({
            "PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
            "PRIMER_NUM_RETURN": self.n_primers,
        })
        if self.primer3_global_args:
            args.update(self.primer3_global_args)
        return args

# -------------------------------------------------------------------------
# 3. Probe 설정 (여기에 필드 추가)
# -------------------------------------------------------------------------
class ProbeKwargs(BaseDesignKwargs):
    # [Override Defaults] Probe는 보통 Primer보다 Tm이 높고 길이가 긴 경향이 있음
    opt_length: int = 25
    min_length: int = 20
    max_length: int = 30
    opt_tm: float = 65.0
    min_tm: float = 60.0
    max_tm: float = 70.0
    min_gc: float = 35.0
    max_gc: float = 65.0
    
    n_probes: int = 10  # Base의 n_primers 대신 사용
    
    # [Tm Diff Logic]
    min_primer_probe_tm_diff: float = 5.0
    max_primer_probe_tm_diff: float = 10.0

    # ✅ [NEW] Python Logic용 Custom Filter Options
    max_probe_poly_g: int = Field(
        4, description="Probe 내 연속된 G의 최대 허용 개수 (예: 4면 GGGG 허용, GGGGG 불가)"
    )
    max_probe_3_end_gc: int = Field(
        4, description="Probe 3' 말단 5bp 내 G/C의 최대 미만 개수 (예: 2면 0, 1개만 허용)"
    )

    def to_global_args(self) -> Dict[str, Any]:
        """Probe 설계를 위한 Global Args 생성 (Prefix: PRIMER_INTERNAL_)"""
        # 1. Base 매핑 (PRIMER_INTERNAL_ 접두어 사용)
        args = self.to_primer3_args(prefix="PRIMER_INTERNAL_")
        
        # 2. Probe 개수 덮어쓰기 (Base의 n_primers 대신 n_probes 사용)
        args["PRIMER_INTERNAL_NUM_RETURN"] = self.n_probes
        
        # 3. Custom Field는 Primer3에 넘기지 않으므로 여기서 제거할 필요는 없지만,
        #    Primer3용 딕셔너리에는 포함되지 않게 함 (위의 args에는 포함 안 됨)
        
        return args

# -------------------------------------------------------------------------
# 4. 기타 설정
# -------------------------------------------------------------------------
class BisulfiteKwargs(BaseModel):
    run: bool = False
    cpg_default: str = "methyl"

class ASPCRKwargs(BaseModel):
    # 기본 설정

    target_allele: Literal["ref", "alt"] = Field("alt", description="Target allele to amplify")
    
    # ----------------------------------------------------------------
    # Artificial Mismatch (AM) 설정
    # ----------------------------------------------------------------
    use_artificial_mismatch: bool = Field(
        False, description="3' 말단 근처에 인위적 불일치(Artificial Mismatch) 도입 여부"
    )
    
    # 위치: n-1 (2nd from 3'), n-2 (3rd from 3')
    # User Notation: n-1 (2nd), n-2 (3rd) -> Python Index: -2, -3
    am_position: Literal["-2", "-3"] = Field(
        "-3", description="인위적 불일치 위치 (-1: 3' 끝에서 두 번째, -2: 세 번째)"
    )
    
    # 강도: Strong (강한 억제), Weak (약한 억제)
    # SNP Mismatch가 약하면(G:T 등) -> Stsrong AM 추천
    # SNP Mismsssatch가 강하면(A:A 등) -> Weak AM 추천
    am_strength: Literal["strong", "weak", "auto"] = Field(
        "auto", description="Mismatch 결합력 억제 강도 (Auto: 3' 말단 염기쌍에 따라 자동 결정)"
    )

    def get_position_index(self) -> int:
        """Python Slicing용 인덱스 반환"""
        return -2 if self.am_position == "n-1" else -3

# -------------------------------------------------------------------------
# 5. 통합 PCR Params
# -------------------------------------------------------------------------
class PCRParams(BaseModel):
    primer_kwargs: PrimerKwargs
    bisulfite_kwargs: BisulfiteKwargs
    as_pcr_kwargs: Optional[ASPCRKwargs] = None
    
    # ✅ Probe는 선택 사항이므로 Optional로 변경
    probe_kwargs: Optional[ProbeKwargs] = None