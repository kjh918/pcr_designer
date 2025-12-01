# app/schemas.py

from typing import Optional, List
from pydantic import BaseModel, Field


class RegionInput(BaseModel):
    chrom: str = Field(..., description="염색체 이름, 예: chr1")
    start: int = Field(..., description="영역 시작 (1-based)")
    end: int = Field(..., description="영역 끝 (1-based, inclusive)")
    name: Optional[str] = Field(None, description="영역 이름/ID (optional)")
    sequence: Optional[str] = Field(
        None,
        description="해당 영역의 염기서열(선택). 없으면 내부에서 FASTA 등으로 가져오도록 구현 가능"
    )


class Primer(BaseModel):
    seq: str = Field(..., description="primer 서열 (5'→3')")
    tm: Optional[float] = Field(None, description="녹는온도(Tm)")
    gc: Optional[float] = Field(None, description="GC 비율(%)")
    start: Optional[int] = Field(None, description="target 내에서 시작 위치")
    end: Optional[int] = Field(None, description="target 내에서 끝 위치")
    strand: Optional[str] = Field(None, description="+ 또는 -")


class PrimerPair(BaseModel):
    forward: Primer
    reverse: Primer
    product_size: Optional[int] = Field(None, description="PCR product size")


class SingleRegionDesignRequest(BaseModel):
    region: RegionInput


class SingleRegionDesignResponse(BaseModel):
    region: RegionInput
    primer_pair: PrimerPair


class MultiRegionDesignResult(BaseModel):
    region: RegionInput
    primer_pair: PrimerPair


class MultiRegionDesignResponse(BaseModel):
    results: List[MultiRegionDesignResult]
