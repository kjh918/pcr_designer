# app/routers/design_manual.py
from __future__ import annotations

from typing import Any, Dict, Optional
import traceback
from datetime import datetime

from fastapi import APIRouter, Request, Form, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from pcr.designers.base import PrimerDesigner, ProbePrimerDesigner
from pcr.pipelines.base import run_pipeline

from pcr.config.runtime import (
    build_pcr_params_from_web,
    build_qc_params_from_web,
    resolve_pcr_params,
    merge_dict,
)

router = APIRouter(prefix="/design", tags=["design"])
templates = Jinja2Templates(directory="app/templates")


def _normalize_seq(s: str) -> str:
    return "".join((s or "").strip().upper().split())


@router.get("/manual", response_class=HTMLResponse, name="manual_page")
async def manual_page(request: Request):
    """
    Manual 입력 폼 페이지 (manual.html)
    """
    return templates.TemplateResponse(
        "manual.html",
        {
            "request": request,
            "assay": "manual",
            "error": None,
        },
    )


# ✅ IMPORTANT: 템플릿이 url_for('design_manual_from_form') 호출하면 이 name이 있어야 함
@router.post("/manual/qpcr", response_class=HTMLResponse, name="design_manual_from_form")
async def design_manual_from_form(
    request: Request,

    # 기본 입력
    name: str = Form(""),
    template_sequence: str = Form(...),
    target_start_index: int = Form(...),
    target_end_index: int = Form(...),

    # qpcr 전용 (probe on/off)
    probe: str = Form("no"),  # yes/no

    # primer 공통 옵션
    min_amplicon_length: Optional[int] = Form(80),
    max_amplicon_length: Optional[int] = Form(120),
    n_primers: Optional[int] = Form(50),

    primer_opt_length: Optional[int] = Form(20),
    primer_min_length: Optional[int] = Form(18),
    primer_max_length: Optional[int] = Form(25),

    primer_opt_tm: Optional[float] = Form(60.0),
    primer_min_tm: Optional[float] = Form(50.0),
    primer_max_tm: Optional[float] = Form(70.0),

    primer_opt_gc: Optional[float] = Form(45.0),
    primer_min_gc: Optional[float] = Form(35.0),
    primer_max_gc: Optional[float] = Form(65.0),

    # probe 옵션
    n_probes: Optional[int] = Form(10),
    probe_opt_length: Optional[int] = Form(25),
    probe_min_length: Optional[int] = Form(20),
    probe_max_length: Optional[int] = Form(30),
    probe_opt_tm: Optional[float] = Form(60.0),
    probe_min_tm: Optional[float] = Form(50.0),
    probe_max_tm: Optional[float] = Form(70.0),
    probe_min_gc: Optional[float] = Form(35.0),
    probe_max_gc: Optional[float] = Form(65.0),
):
    """
    template 기반 qPCR primer(+probe) 디자인 후 design.html에 qpcr 결과 형태로 렌더
    """
    context: Dict[str, Any] = {
        "request": request,
        "mode": "single",
        # ✅ design.html은 assay가 qpcr/methyl/as-pcr 중 하나일 때만 정상 분기함
        "assay": "qpcr",
        "primer_type": "default",
        "reference": "manual",
        "probe": probe,

        "single_result": None,
        "single_total_amplicons": None,
        "single_filtered_amplicons": None,
        "multi_results": None,
        "export_meta": None,
        "error": None,
    }

    try:
        tpl = _normalize_seq(template_sequence)
        if not tpl:
            raise HTTPException(status_code=400, detail="template_sequence is empty")

        ts = int(target_start_index)
        te = int(target_end_index)
        if not (0 <= ts <= te < len(tpl)):
            raise HTTPException(
                status_code=400,
                detail=f"target index out of range: start={ts}, end={te}, len={len(tpl)}",
            )

        # probe on/off 반영
        effective_n_probes = int(n_probes or 0) if (probe == "yes") else 0

        # -----------------------------
        # ✅ 1) web-style config 생성 (qpcr router 흐름과 동일)
        # -----------------------------
        pcr_overrides: Dict[str, Any] = {
            "primer_kwargs": {
                "min_amplicon_length": min_amplicon_length,
                "max_amplicon_length": max_amplicon_length,
                "n_primers": n_primers,

                "opt_length": primer_opt_length,
                "min_length": primer_min_length,
                "max_length": primer_max_length,

                "opt_tm": primer_opt_tm,
                "min_tm": primer_min_tm,
                "max_tm": primer_max_tm,

                "opt_gc": primer_opt_gc,
                "min_gc": primer_min_gc,
                "max_gc": primer_max_gc,
            },
            "probe_kwargs": {
                "n_probes": effective_n_probes,
                "primer3_global_args": {
                    "PRIMER_INTERNAL_OPT_SIZE": probe_opt_length,
                    "PRIMER_INTERNAL_MIN_SIZE": probe_min_length,
                    "PRIMER_INTERNAL_MAX_SIZE": probe_max_length,
                    "PRIMER_INTERNAL_OPT_TM": probe_opt_tm,
                    "PRIMER_INTERNAL_MIN_TM": probe_min_tm,
                    "PRIMER_INTERNAL_MAX_TM": probe_max_tm,
                    "PRIMER_INTERNAL_MIN_GC": probe_min_gc,
                    "PRIMER_INTERNAL_MAX_GC": probe_max_gc,
                },
            },
            "bisulfite": {"run": False},
        }
        pcr_cfg = build_pcr_params_from_web(pcr_overrides)

        qc_overrides: Dict[str, Any] = {
            "MIN_AMP_BP": min_amplicon_length,
            "MAX_AMP_BP": max_amplicon_length,
        }
        qc_params = build_qc_params_from_web(qc_overrides)

        resolved = resolve_pcr_params(
            pcr_cfg=pcr_cfg,
            min_amplicon_length=min_amplicon_length,
            max_amplicon_length=max_amplicon_length,
            n_probes=effective_n_probes,
            n_primers=n_primers,
            bisulfite=False,
        )

        primer3_global_args = merge_dict(base=pcr_cfg.primer_kwargs.primer3_global_args, override=None)
        probe_primer3_global_args = merge_dict(base=pcr_cfg.probe_kwargs.primer3_global_args, override=None)

        # -----------------------------
        # ✅ 2) Designer 선택 (template 기반)
        # -----------------------------
        primer_kwargs = dict(
            template_sequence=tpl,
            reference_template_sequence=tpl,
            target_start_index=ts,
            target_end_index=te,
            min_amplicon_length=int(resolved.min_amplicon_length),
            max_amplicon_length=int(resolved.max_amplicon_length),
            n_primers=int(resolved.n_primers),
            primer3_global_args=primer3_global_args,
        )

        if probe == "yes" and int(resolved.n_probes or 0) > 0:
            designer = ProbePrimerDesigner(
                **primer_kwargs,
                n_probes=int(resolved.n_probes),
                probe_primer3_global_args=probe_primer3_global_args,
            )
        else:
            designer = PrimerDesigner(**primer_kwargs)

        genomic_id = (name or "MANUAL_TEMPLATE").strip()

        # -----------------------------
        # ✅ 3) Pipeline 실행 (qpcr과 동일하게)
        # -----------------------------
        result = run_pipeline(
            genomic_id=genomic_id,
            qc_params=qc_params,
            assay="qpcr",
            designer=designer,
        )

        total_df = result.total_df
        filtered_df = result.filtered_df

        context["single_result"] = {
            "region": {
                "chrom": "TEMPLATE",
                "start": 0,
                "end": len(tpl) - 1,
                "name": genomic_id,
            },
            "total_count": len(total_df),
            "filtered_count": len(filtered_df),
        }
        context["single_total_amplicons"] = total_df.to_dict(orient="records")
        context["single_filtered_amplicons"] = filtered_df.to_dict(orient="records")

        # export meta
        context["export_meta"] = {
            "assay": "manual_qpcr",  # ✅ 구분값
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "reference": "manual_template",
            "region": {
                "chrom": "TEMPLATE",
                "start": 0,
                "end": len(tpl) - 1,
                "name": genomic_id,
            },
            "template": {
                "length": len(tpl),
                "target_start_index": ts,
                "target_end_index": te,
            },
            "probe": (probe == "yes"),
            "pcr_params": {
                "min_amplicon_length": resolved.min_amplicon_length,
                "max_amplicon_length": resolved.max_amplicon_length,
                "n_primers": resolved.n_primers,
                "n_probes": resolved.n_probes,
            },
        }

    except Exception as e:
        traceback.print_exc()
        context["error"] = str(e)

    response = templates.TemplateResponse("design.html", context)
    response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
    response.headers["Pragma"] = "no-cache"
    response.headers["Expires"] = "0"
    return response
