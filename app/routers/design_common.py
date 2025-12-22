# app/routers/design_common.py
from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from fastapi import Form, File, UploadFile, HTTPException

from app.schemas import RegionInput


# -------------------------
# 공통 폼 파라미터 (항상 공통)
# -------------------------
@dataclass
class CommonDesignForm:
    # UI 공통
    mode: str
    primer_type: str
    reference: str

    # region
    chrom: str
    start: Optional[int]
    end: Optional[int]
    name: str

    # primer 공통 옵션
    min_amplicon_length: Optional[int]
    max_amplicon_length: Optional[int]
    n_primers: Optional[int]

    primer_opt_length: Optional[int]
    primer_min_length: Optional[int]
    primer_max_length: Optional[int]

    primer_opt_gc: Optional[float]
    primer_min_gc: Optional[float]
    primer_max_gc: Optional[float]

    # multi
    file: Optional[UploadFile]

    @classmethod
    def as_form(
        cls,
        mode: str = Form("single"),
        primer_type: str = Form("default"),
        reference: str = Form("hg38"),
        chrom: str = Form(""),
        start: int | None = Form(None),
        end: int | None = Form(None),
        name: str = Form(""),
        min_amplicon_length: int | None = Form(None),
        max_amplicon_length: int | None = Form(None),
        n_primers: int | None = Form(None),
        primer_opt_length: int | None = Form(None),
        primer_min_length: int | None = Form(None),
        primer_max_length: int | None = Form(None),
        primer_opt_gc: float | None = Form(None),
        primer_min_gc: float | None = Form(None),
        primer_max_gc: float | None = Form(None),
        file: UploadFile | None = File(None),
    ) -> "CommonDesignForm":
        return cls(
            mode=mode,
            primer_type=primer_type,
            reference=reference,
            chrom=chrom,
            start=start,
            end=end,
            name=name,
            min_amplicon_length=min_amplicon_length,
            max_amplicon_length=max_amplicon_length,
            n_primers=n_primers,
            primer_opt_length=primer_opt_length,
            primer_min_length=primer_min_length,
            primer_max_length=primer_max_length,
            primer_opt_gc=primer_opt_gc,
            primer_min_gc=primer_min_gc,
            primer_max_gc=primer_max_gc,
            file=file,
        )


def build_common_kwargs(f: CommonDesignForm) -> Dict[str, Any]:
    """pipeline 함수에 공통으로 넘길 kwargs"""
    return dict(
        min_amplicon_length=f.min_amplicon_length,
        max_amplicon_length=f.max_amplicon_length,
        n_primers=f.n_primers,
        primer_opt_length=f.primer_opt_length,
        primer_min_length=f.primer_min_length,
        primer_max_length=f.primer_max_length,
        primer_opt_gc=f.primer_opt_gc,
        primer_min_gc=f.primer_min_gc,
        primer_max_gc=f.primer_max_gc,
    )


async def parse_regions_from_form(f: CommonDesignForm) -> List[RegionInput]:
    """
    mode=single -> 단일 RegionInput 1개 반환
    mode=multi  -> 엑셀에서 여러 RegionInput 반환
    """
    if f.file is not None and getattr(f.file, "filename", ""):
        contents = await f.file.read()
        df = pd.read_excel(BytesIO(contents))

        required_cols = ["chrom", "start", "end", "name"]
        for col in required_cols:
            if col not in df.columns:
                raise HTTPException(status_code=400, detail=f"필수 컬럼이 없습니다: {col}")

        regions: List[RegionInput] = []
        for _, row in df.iterrows():
            chrom_val = str(row["chrom"])
            start_val = int(row["start"])
            end_val = int(row["end"])
            name_val = (
                str(row["name"])
                if "name" in df.columns and not pd.isna(row["name"])
                else None
            )
            regions.append(
                RegionInput(
                    chrom=chrom_val,
                    start=start_val,
                    end=end_val,
                    name=name_val,
                    sequence="",
                )
            )
        return regions

    # ✅ file이 없으면 mode로 판단
    if f.mode == "single":
        if not f.chrom or f.start is None or f.end is None:
            raise HTTPException(status_code=400, detail="chrom / start / end 가 필요합니다.")
        return [
            RegionInput(
                chrom=f.chrom,
                start=f.start,
                end=f.end,
                name=(f.name or None),
                sequence="",
            )
        ]

    if f.mode == "multi":
        if f.file is None or not f.file.filename:
            raise HTTPException(status_code=400, detail="Multiple 모드에서는 Excel 파일이 필요합니다.")

        contents = await f.file.read()
        df = pd.read_excel(BytesIO(contents))

        required_cols = ["chrom", "start", "end", "name"]
        for col in required_cols:
            if col not in df.columns:
                raise HTTPException(status_code=400, detail=f"필수 컬럼이 없습니다: {col}")

        regions: List[RegionInput] = []
        for _, row in df.iterrows():
            chrom_val = str(row["chrom"])
            start_val = int(row["start"])
            end_val = int(row["end"])
            name_val = (
                str(row["name"])
                if "name" in df.columns and not pd.isna(row["name"])
                else None
            )
            regions.append(
                RegionInput(
                    chrom=chrom_val,
                    start=start_val,
                    end=end_val,
                    name=name_val,
                    sequence="",
                )
            )
        return regions

    raise HTTPException(status_code=400, detail=f"알 수 없는 mode: {f.mode}")


def init_context(request, f: CommonDesignForm, *, assay: str) -> Dict[str, Any]:
    """템플릿 공통 context 생성"""
    return {
        "request": request,
        "assay": assay,  # qpcr / methyl / as-pcr
        "mode": f.mode,
        "primer_type": f.primer_type,
        "reference": f.reference,
        "single_result": None,
        "single_total_amplicons": None,
        "single_filtered_amplicons": None,
        "multi_results": None,
        "error": None,
    }
