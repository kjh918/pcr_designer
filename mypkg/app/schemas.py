from typing import Optional, List
from pydantic import BaseModel, Field

class RegionInput(BaseModel):
    chrom: str = Field(..., description="염색체 이름, 예: chr1")
    start: int = Field(..., description="영역 시작 (1-based)")
    end: int = Field(..., description="영역 끝 (1-based, inclusive)")


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
    probe: Optional[Primer] = Field(
        None,
        description="qPCR probe (있는 경우)",
    )
    product_size: Optional[int] = Field(None, description="PCR product size")


# ===============================
#   Primer / Probe 디자인 옵션
#   (웹에서 override 값 받기용)
# ===============================

class PrimerDesignOptions(BaseModel):
    # config.pcr_params.primer_kwargs 와 1:1 매핑
    min_amplicon_length: Optional[int] = Field(
        None, description="최소 amplicon 길이 (없으면 config 기본값)"
    )
    max_amplicon_length: Optional[int] = Field(
        None, description="최대 amplicon 길이 (없으면 config 기본값)"
    )
    n_primers: Optional[int] = Field(
        None, description="설계할 primer 후보 개수 (없으면 config 기본값)"
    )

    opt_length: Optional[int] = Field(
        None, description="Primer 최적 길이 (없으면 config 기본값)"
    )
    min_length: Optional[int] = Field(
        None, description="Primer 최소 길이 (없으면 config 기본값)"
    )
    max_length: Optional[int] = Field(
        None, description="Primer 최대 길이 (없으면 config 기본값)"
    )

    opt_gc: Optional[float] = Field(
        None, description="Primer 최적 GC% (없으면 config 기본값)"
    )
    min_gc: Optional[float] = Field(
        None, description="Primer 최소 GC% (없으면 config 기본값)"
    )
    max_gc: Optional[float] = Field(
        None, description="Primer 최대 GC% (없으면 config 기본값)"
    )

    opt_tm: Optional[float] = Field(
        None, description="Primer 최적 Tm (없으면 config 기본값)"
    )
    min_tm: Optional[float] = Field(
        None, description="Primer 최소 Tm (없으면 config 기본값)"
    )
    max_tm: Optional[float] = Field(
        None, description="Primer 최대 Tm (없으면 config 기본값)"
    )

    # 필요하면 opt_tm / min_tm / max_tm 도 여기서 확장 가능


class ProbeDesignOptions(BaseModel):
    # 웹에서 probe = yes/no 토글 → 여기서는 bool로
    use_probe: Optional[bool] = Field(
        None,
        description="Probe 사용 여부 (웹에서 Probe = Yes/No 토글). None이면 config/biz 로직에서 결정",
    )

    # config.pcr_params.probe_kwargs 와 1:1 매핑
    n_probes: Optional[int] = Field(
        None, description="설계할 probe 개수 (없으면 config 기본값 또는 use_probe=false 시 0)"
    )
    min_primer_probe_tm_diff: Optional[float] = Field(
        None, description="Primer–Probe 최소 Tm 차이"
    )
    max_primer_probe_tm_diff: Optional[float] = Field(
        None, description="Primer–Probe 최대 Tm 차이"
    )

    opt_length: Optional[int] = Field(
        None, description="Probe 최적 길이 (없으면 config 기본값)"
    )
    min_length: Optional[int] = Field(
        None, description="Probe 최소 길이 (없으면 config 기본값)"
    )
    max_length: Optional[int] = Field(
        None, description="Probe 최대 길이 (없으면 config 기본값)"
    )

    opt_tm: Optional[float] = Field(
        None, description="Probe 최적 Tm (없으면 config 기본값)"
    )
    min_tm: Optional[float] = Field(
        None, description="Probe 최소 Tm (없으면 config 기본값)"
    )
    max_tm: Optional[float] = Field(
        None, description="Probe 최대 Tm (없으면 config 기본값)"
    )

    opt_gc: Optional[float] = Field(
        None, description="Probe 최적 GC% (없으면 config 기본값)"
    )
    min_gc: Optional[float] = Field(
        None, description="Probe 최소 GC% (없으면 config 기본값)"
    )
    max_gc: Optional[float] = Field(
        None, description="Probe 최대 GC% (없으면 config 기본값)"
    )


# ===============================
#   요청 / 응답 스키마
# ===============================

class SingleRegionDesignRequest(BaseModel):
    region: RegionInput
    reference: Optional[str] = Field(
        "hg19",
        description="reference 이름 (예: hg19, hg38)",
    )
    primer_options: Optional[PrimerDesignOptions] = Field(
        None,
        description="Primer 설계 옵션 (없으면 config 기본값 사용)",
    )
    probe_options: Optional[ProbeDesignOptions] = Field(
        None,
        description="Probe 설계 옵션 (없으면 config 기본값/Probe 미사용)",
    )


class SingleRegionDesignResponse(BaseModel):
    region: RegionInput
    primer_pair: PrimerPair


class MultiRegionDesignResult(BaseModel):
    region: RegionInput
    primer_pair: PrimerPair


class MultiRegionDesignResponse(BaseModel):
    results: List[MultiRegionDesignResult]


class MultiRegionDesignRequest(BaseModel):
    regions: List[RegionInput]
    reference: Optional[str] = Field(
        "hg19",
        description="reference 이름 (예: hg19, hg38)",
    )
    primer_options: Optional[PrimerDesignOptions] = Field(
        None,
        description="Primer 설계 옵션 (없으면 config 기본값 사용)",
    )
    probe_options: Optional[ProbeDesignOptions] = Field(
        None,
        description="Probe 설계 옵션 (없으면 config 기본값/Probe 미사용)",
    )
