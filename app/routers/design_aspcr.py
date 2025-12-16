# app/routers/design_aspcr.py
from __future__ import annotations

from typing import Any, Dict, List
import traceback

from fastapi import APIRouter, Request, Depends, Form, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from app.routers.design_common import (
    CommonDesignForm,
    build_common_kwargs,
    parse_regions_from_form,
    init_context,
)

# as-pcr pipeline 함수가 있다고 가정
# from pcr.pipelines.variant_probe import design_aspcr_for_region
from pcr.pipelines.as_pcr import run_as_pcr_pipeline  # 임시


router = APIRouter(prefix="/design", tags=["design"])
templates = Jinja2Templates(directory="app/templates")


@router.post("/as-pcr", response_class=HTMLResponse)
async def design_aspcr_from_form(
    request: Request,
    f: CommonDesignForm = Depends(CommonDesignForm.as_form),

    # AS-PCR 전용: variant 정보(예시)
    variant_index: int | None = Form(None),  # template 기준 index or genomic pos 등 네 정책에 맞춰
    ref_allele: str = Form(""),
    alt_allele: str = Form(""),
):
    context = init_context(request, f, assay="as-pcr")

    try:
        if variant_index is None or not ref_allele or not alt_allele:
            raise HTTPException(status_code=400, detail="AS-PCR은 variant_index/ref_allele/alt_allele가 필요합니다.")

        regions = await parse_regions_from_form(f)

        kwargs = build_common_kwargs(f)
        kwargs.update(
            dict(
                bisulfite=False,
                # as-pcr 전용 인자들 (pipeline에 맞게 이름 정리)
                variant_index=variant_index,
                ref_allele=ref_allele,
                alt_allele=alt_allele,
            )
        )

        # 임시로 qpcr 호출. 실제로는 variant pipeline 함수로 대체하면 됨.
        if f.mode == "single":
            region = regions[0]
            total_df, filtered_df = design_qpcr_for_region(
                region=region, reference_name=f.reference, **kwargs
            )
            context["single_result"] = {
                "region": region,
                "total_count": len(total_df),
                "filtered_count": len(filtered_df),
            }
            context["single_total_amplicons"] = total_df.to_dict(orient="records")
            context["single_filtered_amplicons"] = filtered_df.to_dict(orient="records")
        else:
            multi_results: List[Dict[str, Any]] = []
            for region in regions:
                total_df, filtered_df = design_qpcr_for_region(
                    region=region, reference_name=f.reference, **kwargs
                )
                multi_results.append(
                    dict(
                        region=region,
                        total_count=len(total_df),
                        filtered_count=len(filtered_df),
                        total_amplicons=total_df.to_dict(orient="records"),
                        filtered_amplicons=filtered_df.to_dict(orient="records"),
                    )
                )
            context["multi_results"] = multi_results

    except Exception as e:
        traceback.print_exc()
        context["error"] = str(e)

    response = templates.TemplateResponse("design.html", context)
    response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
    response.headers["Pragma"] = "no-cache"
    response.headers["Expires"] = "0"
    return response
