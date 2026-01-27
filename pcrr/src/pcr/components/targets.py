from dataclasses import dataclass
from typing import Literal, Union

# ----------------------------------------------------------------
# 1. Variant Target (SNP, INDEL)
# ----------------------------------------------------------------
@dataclass
class Variant:
    chrom: str
    start: int        # Genomic Start (1-based)
    end: int          # Genomic End (1-based, inclusive)
    id: str           # rsID etc.
    ref: str
    alt: str
    
    # 템플릿 서열 내에서의 좌표 (0-based, Python Slicing style)
    # SNP라면: index_end = index_start + 1
    # Deletion라면: ref 길이만큼 차지
    # Insertion라면: alt 길이만큼 차지 (Alt Template 기준)
    index_start: int  
    index_end: int    
    
    # Target Type
    target_allele: Literal["ref", "alt"] = "alt"

    def get_expected_sequence(self) -> str:
        """Probe 서열 상에서 기대되는 전체 서열 (String)"""
        # 템플릿이 이미 해당 Allele로 치환되어 있다고 가정
        return self.ref if self.target_allele == "ref" else self.alt


# ----------------------------------------------------------------
# 2. CpG Target (Methylation)
# ----------------------------------------------------------------
@dataclass
class CpG:
    chrom: str
    pos: int          # Genomic Coordinate (C position)
    index: int        # 0-based index of 'C'
    
    methylation_status: Literal["methylated", "unmethylated"] = "methylated"

    # Variant와 인터페이스 통일을 위해 property 사용
    @property
    def index_start(self) -> int:
        return self.index

    @property
    def index_end(self) -> int:
        return self.index + 1  # 1 base length

    def get_expected_sequence(self) -> str:
        """
        Bisulfite 처리된 템플릿 기준 기대 염기 반환.
        """
        if self.methylation_status == "methylated":
            return "C"  # Methylated C -> C
        else:
            return "T"  # Unmethylated C -> T


# ----------------------------------------------------------------
# 3. Union Type
# ----------------------------------------------------------------
TargetType = Union[Variant, CpG]