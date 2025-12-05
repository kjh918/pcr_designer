# app/routers/design.py

from __future__ import annotations

from io import BytesIO
from typing import Any, Dict, List, Optional

import json, openpyxl
import pandas as pd
from fastapi import APIRouter, Request, Form, File, UploadFile, HTTPException
from fastapi.responses import HTMLResponse, StreamingResponse
from fastapi.templating import Jinja2Templates

from app.schemas import RegionInput
from primer.core import design_qpcr_for_region  # qPCR wrapper

router = APIRouter(
    prefix="/design",      # ★ /design/xxx 로 묶인다
    tags=["design"],       # Swagger 문서에서 그룹 이름
)

templates = Jinja2Templates(directory="app/templates")


def build_design_kwargs_from_form(
    *,
    probe: str,
    min_amplicon_length: Optional[int],
    max_amplicon_length: Optional[int],
    n_primers: Optional[int],
    n_probes: Optional[int],
    primer_opt_length: Optional[int],
    primer_min_length: Optional[int],
    primer_max_length: Optional[int],
    primer_opt_gc: Optional[float],
    primer_min_gc: Optional[float],
    primer_max_gc: Optional[float],
    min_primer_probe_tm_diff: Optional[float],
    max_primer_probe_tm_diff: Optional[float],
    probe_opt_length: Optional[int],
    probe_min_length: Optional[int],
    probe_max_length: Optional[int],
    probe_opt_tm: Optional[float],
    probe_min_tm: Optional[float],
    probe_max_tm: Optional[float],
    probe_opt_gc: Optional[float],
    probe_min_gc: Optional[float],
    probe_max_gc: Optional[float],
    bisulfite: Optional[bool] = None,
) -> Dict[str, Any]:
    # probe yes/no → n_probes 정리
    effective_n_probes = n_probes if probe == "yes" else 0

    return dict(
        min_amplicon_length=min_amplicon_length,
        max_amplicon_length=max_amplicon_length,
        n_primers=n_primers,
        n_probes=effective_n_probes,
        bisulfite=bisulfite,
        primer_opt_length=primer_opt_length,
        primer_min_length=primer_min_length,
        primer_max_length=primer_max_length,
        primer_opt_gc=primer_opt_gc,
        primer_min_gc=primer_min_gc,
        primer_max_gc=primer_max_gc,
        min_primer_probe_tm_diff=min_primer_probe_tm_diff,
        max_primer_probe_tm_diff=max_primer_probe_tm_diff,
        probe_opt_length=probe_opt_length,
        probe_min_length=probe_min_length,
        probe_max_length=probe_max_length,
        probe_opt_tm=probe_opt_tm,
        probe_min_tm=probe_min_tm,
        probe_max_tm=probe_max_tm,
        probe_opt_gc=probe_opt_gc,
        probe_min_gc=probe_min_gc,
        probe_max_gc=probe_max_gc,
    )


@router.post("/form", response_class=HTMLResponse)
async def design_from_form(
    request: Request,
    mode: str = Form("single"),
    primer_type: str = Form("default"),
    reference: str = Form("hg19"),
    probe: str = Form("no"),
    methylation: str = Form("no"),        # yes / no → bisulfite 플래그

    chrom: str = Form(""),
    start: int | None = Form(None),
    end: int | None = Form(None),
    name: str = Form(""),

    # Primer 설정
    min_amplicon_length: int | None = Form(None),
    max_amplicon_length: int | None = Form(None),
    n_primers: int | None = Form(None),
    primer_opt_length: int | None = Form(None),
    primer_min_length: int | None = Form(None),
    primer_max_length: int | None = Form(None),
    primer_opt_gc: float | None = Form(None),
    primer_min_gc: float | None = Form(None),
    primer_max_gc: float | None = Form(None),

    # Probe 설정
    n_probes: int | None = Form(None),
    min_primer_probe_tm_diff: float | None = Form(None),
    max_primer_probe_tm_diff: float | None = Form(None),
    probe_opt_length: int | None = Form(None),
    probe_min_length: int | None = Form(None),
    probe_max_length: int | None = Form(None),
    probe_opt_tm: float | None = Form(None),
    probe_min_tm: float | None = Form(None),
    probe_max_tm: float | None = Form(None),
    probe_opt_gc: float | None = Form(None),
    probe_min_gc: float | None = Form(None),
    probe_max_gc: float | None = Form(None),

    file: UploadFile | None = File(None),
):
    context: Dict[str, Any] = {
        "request": request,
        "mode": mode,
        "primer_type": primer_type,
        "reference": reference,
        "probe": probe,
        "methylation": methylation,
        "single_result": None,
        "single_total_amplicons": None,
        "single_filtered_amplicons": None,
        "multi_results": None,
        "error": None,
    }

    try:
        ## TODO 
        bisulfite_flag = methylation == "yes"

        design_kwargs = build_design_kwargs_from_form(
            probe=probe,
            min_amplicon_length=min_amplicon_length,
            max_amplicon_length=max_amplicon_length,
            n_primers=n_primers,
            n_probes=n_probes,
            primer_opt_length=primer_opt_length,
            primer_min_length=primer_min_length,
            primer_max_length=primer_max_length,
            primer_opt_gc=primer_opt_gc,
            primer_min_gc=primer_min_gc,
            primer_max_gc=primer_max_gc,
            min_primer_probe_tm_diff=min_primer_probe_tm_diff,
            max_primer_probe_tm_diff=max_primer_probe_tm_diff,
            probe_opt_length=probe_opt_length,
            probe_min_length=probe_min_length,
            probe_max_length=probe_max_length,
            probe_opt_tm=probe_opt_tm,
            probe_min_tm=probe_min_tm,
            probe_max_tm=probe_max_tm,
            probe_opt_gc=probe_opt_gc,
            probe_min_gc=probe_min_gc,
            probe_max_gc=probe_max_gc,
            bisulfite=bisulfite_flag,
        )

        # 1) Single 모드
        if mode == "single":
            if not chrom or start is None or end is None:
                raise HTTPException(
                    status_code=400,
                    detail="qPCR 설계에는 chrom / start / end 가 모두 필요합니다.",
                )

            region = RegionInput(
                chrom=chrom,
                start=start,
                end=end,
                name=name or None,
                sequence="",
            )

            total_df, filtered_df = design_qpcr_for_region(
                region=region,
                reference_name=reference,
                **design_kwargs,
            )

            total_records = total_df.to_dict(orient="records")
            filtered_records = filtered_df.to_dict(orient="records")

            context["single_result"] = {
                "region": region,
                "total_count": len(total_records),
                "filtered_count": len(filtered_records),
            }
            context["single_total_amplicons"] = total_records
            context["single_filtered_amplicons"] = filtered_records
            print(filtered_records)

        # 2) Multi 모드
        elif mode == "multi":
            if file is None or not file.filename:
                raise HTTPException(
                    status_code=400,
                    detail="Multiple 모드에서는 Excel 파일이 필요합니다.",
                )

            contents = await file.read()
            df = pd.read_excel(BytesIO(contents))

            required_cols = ["chrom", "start", "end"]
            for col in required_cols:
                if col not in df.columns:
                    raise HTTPException(
                        status_code=400,
                        detail=f"필수 컬럼이 없습니다: {col}",
                    )

            multi_results: List[Dict[str, Any]] = []

            for _, row in df.iterrows():
                chrom_val = str(row["chrom"])
                start_val = int(row["start"])
                end_val = int(row["end"])
                name_val = (
                    str(row["name"])
                    if "name" in df.columns and not pd.isna(row["name"])
                    else None
                )

                region = RegionInput(
                    chrom=chrom_val,
                    start=start_val,
                    end=end_val,
                    name=name_val,
                    sequence="",
                )

                total_df, filtered_df = design_qpcr_for_region(
                    region=region,
                    reference_name=reference,
                    **design_kwargs,
                )

                multi_results.append(
                    {
                        "region": region,
                        "total_count": len(total_df),
                        "filtered_count": len(filtered_df),
                        "total_amplicons": total_df.to_dict(orient="records"),
                        "filtered_amplicons": filtered_df.to_dict(orient="records"),
                    }
                )

            context["multi_results"] = multi_results

        else:
            raise HTTPException(
                status_code=400,
                detail=f"알 수 없는 mode: {mode}",
            )

    except HTTPException as e:
        context["error"] = e.detail
    except ValueError as e:
        context["error"] = f"Reference/설계 에러: {e}"
    except Exception as e:
        context["error"] = str(e)

    return templates.TemplateResponse("index.html", context)


@router.post("/export")
async def export_amplicons_to_excel(
    kind: str = Form(...),      # "single_total" / "single_filtered" / "multi_..." 등
    data_json: str = Form(...), # 템플릿에서 넘겨주는 JSON 문자열
):
    """
    템플릿에서 넘겨준 amplicon 리스트(JSON)를 엑셀로 변환하여 다운로드.
    """
    try:
        rows = json.loads(data_json)
        if not isinstance(rows, list):
            raise ValueError("data_json must be a list of dicts")

        df = pd.DataFrame(rows)
        print(df)
        buffer = BytesIO()
        with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
            sheet_name = kind[:31] if kind else "Sheet1"  # 엑셀 시트 이름 최대 31자
            df.to_excel(writer, index=False, sheet_name=sheet_name)

        buffer.seek(0)

        filename = f"{kind}_amplicons.xlsx"

        return StreamingResponse(
            buffer,
            media_type=(
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            ),
            headers={
                "Content-Disposition": f'attachment; filename="{filename}"'
            },
        )
    except Exception as e:
        # 프론트용 메시지로 쓰려면 TemplateResponse로 돌려도 되는데,
        # 여기서는 단순 HTTP 400 에러로 처리
        raise HTTPException(status_code=400, detail=f"엑셀 export 실패: {e}")
