# app/routers/design_qpcr.py
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

from pcr.seq.fetch import GenomicRegion
from pcr.pipelines.qpcr import run_qpcr

from pcr.config.runtime import settings as get_settings, get_fasta_handle, resolve_pcr_params, merge_dict


router = APIRouter(prefix="/design", tags=["design"])
templates = Jinja2Templates(directory="app/templates")


@router.post("/qpcr", response_class=HTMLResponse)
async def design_qpcr_from_form(
    request: Request,
    f: CommonDesignForm = Depends(CommonDesignForm.as_form),

    # qpcr 전용
    probe: str = Form("no"),  # yes/no

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
):
    context = init_context(request, f, assay="qpcr")

    try:
        s = get_settings()
        pcr_cfg = s.pcr_params
        qc_cfg = s.qc_params

        # FASTA handle (pysam cached)
        fasta = get_fasta_handle(f.reference)

        # region(s) 파싱
        regions = await parse_regions_from_form(f)

        # 공통 kwargs (라우터 수준)
        common_kwargs = build_common_kwargs(f)

        # probe on/off 반영
        effective_n_probes = (n_probes if n_probes is not None else None)
        if probe != "yes":
            effective_n_probes = 0

        # settings 기본값 + override 확정
        resolved = resolve_pcr_params(
            min_amplicon_length=common_kwargs.get("min_amplicon_length"),
            max_amplicon_length=common_kwargs.get("max_amplicon_length"),
            n_probes=effective_n_probes,
            n_primers=common_kwargs.get("n_primers"),
            bisulfite=False,
        )

        # primer3 args는 config 기본 + (원하면 나중에 form에서 override) merge
        primer3_global_args = merge_dict(
            base=pcr_cfg.primer_kwargs.primer3_global_args,
            override=None,
        )
        probe_primer3_global_args = merge_dict(
            base=pcr_cfg.probe_kwargs.primer3_global_args,
            override=None,
        )

        # 결과 조립
        if f.mode == "single":
            region = regions[0]
            gr = GenomicRegion(chrom=region.chrom, start=region.start, end=region.end)

            result = run_qpcr(
                region=gr,
                fasta=fasta,
                pcr_cfg=pcr_cfg,
                qc_cfg=qc_cfg,
                min_amplicon_length=resolved.min_amplicon_length,
                max_amplicon_length=resolved.max_amplicon_length,
                n_probes=resolved.n_probes,
                n_primers=resolved.n_primers,
                primer3_global_args=primer3_global_args,
                probe_primer3_global_args=probe_primer3_global_args,
                # TODO: 아래 primer/probe override들 연결하려면 run_qpcr 시그니처 확장
                # primer_opt_length=f.primer_opt_length, ...
            )

            # PipelineResult → 템플릿용 dict
            total_df = result.total_df
            filtered_df = result.filtered_df

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
                gr = GenomicRegion(chrom=region.chrom, start=region.start, end=region.end)

                result = run_qpcr(
                    region=gr,
                    fasta=fasta,
                    pcr_cfg=pcr_cfg,
                    qc_cfg=qc_cfg,
                    min_amplicon_length=resolved.min_amplicon_length,
                    max_amplicon_length=resolved.max_amplicon_length,
                    n_probes=resolved.n_probes,
                    n_primers=resolved.n_primers,
                    primer3_global_args=primer3_global_args,
                    probe_primer3_global_args=probe_primer3_global_args,
                )

                total_df = result.total_df
                filtered_df = result.filtered_df

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
