# app/routers/qc.py

from __future__ import annotations

from fastapi import APIRouter, Request, Form
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from primer.qc import (
    QCThresholds,
    evaluate_amplicons,
    make_amplicon_for_qc,
    blast_qc_for_primer_pair,
)
from config.settings import settings

templates = Jinja2Templates(directory="app/templates")

router = APIRouter(
    prefix="/qc",
    tags=["qc"],
)


@router.get("/", response_class=HTMLResponse)
async def qc_page(
    request: Request,
    error: str | None = None,
):
    """
    QC Only 화면: 처음에는 결과 없이 폼만 보여줌
    """
    qc_result = None

    return templates.TemplateResponse(
        "qc.html",   # ✅ index.html → qc.html 로 변경
        {
            "request": request,
            "error": error,
            "qc_result": qc_result,

            # 폼 기본값들 (템플릿에서 value로 쓸 수 있게)
            "forward": "",
            "reverse": "",
            "probe": "",
            "template_sequence": "",
            "reference": "hg38",
        },
    )


@router.post("/only", response_class=HTMLResponse)
async def qc_only_run(
    request: Request,
    forward: str = Form(...),
    reverse: str = Form(...),
    probe: str | None = Form(None),
    template_sequence: str | None = Form(None),
    reference: str = Form("hg38"),
):
    """
    Primer / Probe 서열만 넣어서
    - primer3 기반 Thermo QC (hairpin / homodimer / heterodimer)
    - BLAST 기반 off-target / self-amplicon QC
    를 한 번에 계산하는 엔드포인트.
    """
    error: str | None = None
    qc_result = None

    try:
        # 1) Thermo 기반 QC
        amp = make_amplicon_for_qc(
            forward_seq=forward,
            reverse_seq=reverse,
            probe_seq=probe,
            template_seq=template_sequence,
        )

        qc_thresholds = QCThresholds()  # 기본값은 settings.qc_params에서 온 값 사용
        total_rows, filtered_rows = evaluate_amplicons(
            genomic_id="QC_ONLY",
            amplicons=[amp],
            qc_thresholds=qc_thresholds,
        )

        # 하나만 평가했으므로 첫번째 row 사용
        thermo_row = total_rows[0]

        # 2) BLAST 기반 QC
        ref_cfg = settings.references[reference]
        blast_db = str(ref_cfg.blast)

        blast_qc = blast_qc_for_primer_pair(
            f_name="FORWARD",
            f_seq=forward,
            r_name="REVERSE",
            r_seq=reverse,
            db=blast_db,
        )

        # 3) Thermo + BLAST 결과 합치기
        qc_result = {
            **thermo_row,
            "BLAST_F_HITS": blast_qc["f_hits"],
            "BLAST_R_HITS": blast_qc["r_hits"],
            "BLAST_NEARBY_AMP_COUNT": blast_qc["nearby_count"],
            "BLAST_MIN_AMP_SIZE": blast_qc["min_amplicon_size"],
            "BLAST_AMP_DETAIL": ";".join(blast_qc["amplicon_details"])
            if blast_qc["amplicon_details"]
            else "",
            "QC_BLAST_HIT": blast_qc["qc_blast_hit"],
            "QC_BLAST_AMP": blast_qc["qc_blast_amplicon"],
        }

    except Exception as e:
        error = f"QC 수행 중 오류: {e}"

    return templates.TemplateResponse(
        "qc.html",   # ✅ 여기서도 qc.html 렌더
        {
            "request": request,
            "error": error,
            "qc_result": qc_result,

            # 폼에 다시 채워줄 값들 (에러나도 입력 유지)
            "forward": forward,
            "reverse": reverse,
            "probe": probe or "",
            "template_sequence": template_sequence or "",
            "reference": reference,
        },
    )