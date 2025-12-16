# app/routers/design_methyl.py
from __future__ import annotations

from typing import Any, Dict, List
import traceback

from fastapi import APIRouter, Request, Depends, Form
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from app.routers.design_common import (
    CommonDesignForm,
    build_common_kwargs,
    parse_regions_from_form,
    init_context,
)

# methyl용 pipeline 함수가 따로 있다고 가정 (없으면 qpcr에 bisulfite=True로 연결해도 됨)
from pcr.pipelines.qpcr import run_qpcr


router = APIRouter(prefix="/design", tags=["design"])
templates = Jinja2Templates(directory="app/templates")


@router.post("/methyl", response_class=HTMLResponse)
async def design_methyl_from_form(
    request: Request,
    f: CommonDesignForm = Depends(CommonDesignForm.as_form),

    # methyl 전용 (필요한 만큼만)
    cpg_default: str = Form("methyl"),  # parameter.yaml 기본값과 맞추면 좋음
):
    context = init_context(request, f, assay="methyl")

    try:
        regions = await parse_regions_from_form(f)

        kwargs = build_common_kwargs(f)
        kwargs.update(
            dict(
                bisulfite=True,
                cpg_default=cpg_default,  # pipeline/designer가 받는다면 전달
                # methyl은 보통 probe 옵션이 달라질 수 있으니 여기선 qpcr probe 인자 아예 안 받음
                n_probes=0,
            )
        )

        if f.mode == "single":
            region = regions[0]
            total_df, filtered_df = run_qpcr(
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
                total_df, filtered_df = run_qpcr(
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
