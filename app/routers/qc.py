# app/routers/qc.py
# app/routers/qc.py
from __future__ import annotations

from fastapi import APIRouter, Request, Form
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from pcr.config.settings import settings
from pcr.config.runtime import build_qc_params_from_web

from pcr.qc.evaluate import run_thermo_qc, run_blast_qc
from pcr.qc.factories import make_amplicon_for_qc

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
    QC Only 화면 (초기 상태)
    """
    return templates.TemplateResponse(
        "qc.html",
        {
            "request": request,
            "error": error,
            "qc_result": None,

            # form 기본값
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

    # --- primer / probe ---
    forward: str = Form(...),
    reverse: str = Form(...),
    probe: str | None = Form(None),
    template_sequence: str | None = Form(None),

    # --- reference ---
    reference: str = Form("hg38"),

    # =========================
    # QC override (대문자 기준)
    # =========================
    PRIMER_MAX_DIFF_TM: float | None = Form(None),
    PRIMER_MIN_DIFF_TM: float | None = Form(None),
    PROBE_MAX_DIFF_TM: float | None = Form(None),
    PROBE_MIN_DIFF_TM: float | None = Form(None),

    HAIRPIN_MIN_DG: float | None = Form(None),
    HOMODIMER_MIN_DG: float | None = Form(None),
    HETERODIMER_MIN_DG: float | None = Form(None),

    BLAST_IDENTITY_THRESHOLD: float | None = Form(None),
    BLAST_LENGTH_THRESHOLD: int | None = Form(None),
    BLAST_MAX_ALIGNMENTS: int | None = Form(None),
    MIN_AMP_BP: int | None = Form(None),
    MAX_AMP_BP: int | None = Form(None),
):
    """
    Primer / Probe sequence만으로
    - Thermo QC (primer3)
    - BLAST QC (off-target / self-amplicon)
    실행
    """
    error: str | None = None
    qc_result = None

    try:
        # -------------------------------------------------
        # 1) QCParams 생성 (settings + web override)
        # -------------------------------------------------
        qc_params = build_qc_params_from_web({
            "PRIMER_MAX_DIFF_TM": PRIMER_MAX_DIFF_TM,
            "PRIMER_MIN_DIFF_TM": PRIMER_MIN_DIFF_TM,
            "PROBE_MAX_DIFF_TM": PROBE_MAX_DIFF_TM,
            "PROBE_MIN_DIFF_TM": PROBE_MIN_DIFF_TM,
            "HAIRPIN_MIN_DG": HAIRPIN_MIN_DG,
            "HOMODIMER_MIN_DG": HOMODIMER_MIN_DG,
            "HETERODIMER_MIN_DG": HETERODIMER_MIN_DG,
            "BLAST_IDENTITY_THRESHOLD": BLAST_IDENTITY_THRESHOLD,
            "BLAST_LENGTH_THRESHOLD": BLAST_LENGTH_THRESHOLD,
            "BLAST_MAX_ALIGNMENTS": BLAST_MAX_ALIGNMENTS,
            "MIN_AMP_BP": MIN_AMP_BP,
            "MAX_AMP_BP": MAX_AMP_BP,
        })

        # -------------------------------------------------
        # 2) Thermo QC
        # -------------------------------------------------
        amp = make_amplicon_for_qc(
            forward_seq=forward,
            reverse_seq=reverse,
            probe_seq=probe,
            template_seq=template_sequence,
        )

        total_rows, _filtered_rows = run_thermo_qc(
            genomic_id="QC_ONLY",
            amplicons=[amp],
            qc_params=qc_params,
        )

        thermo_row = total_rows[0] if total_rows else {}

        # -------------------------------------------------
        # 3) BLAST QC
        # -------------------------------------------------
        if reference not in settings.references:
            raise ValueError(f"Unknown reference: {reference}")

        blast_db = str(settings.references[reference].blast)

        blast_qc = run_blast_qc(
            f_name="FORWARD",
            f_seq=forward,
            r_name="REVERSE",
            r_seq=reverse,
            db=blast_db,
            qc_params=qc_params,
            probe_name="PROBE" if probe else None,
            probe_seq=probe,
        )

        # -------------------------------------------------
        # 4) 결과 병합
        # -------------------------------------------------
        qc_result = {
            **thermo_row,
            "BLAST_F_HITS": blast_qc.get("f_hits"),
            "BLAST_R_HITS": blast_qc.get("r_hits"),
            "BLAST_NEARBY_AMP_COUNT": blast_qc.get("nearby_count"),
            "BLAST_MIN_AMP_SIZE": blast_qc.get("min_amplicon_size"),
            "BLAST_AMP_DETAIL": ";".join(blast_qc.get("amplicon_details") or []),
            "QC_BLAST_HIT": blast_qc.get("qc_blast_hit"),
            "QC_BLAST_AMP": blast_qc.get("qc_blast_amplicon"),
            "QC_PROBE_IN_AMP": blast_qc.get("qc_probe_in_amplicon", "-"),
        }

    except Exception as e:
        error = f"QC 수행 중 오류: {e}"

    return templates.TemplateResponse(
        "qc.html",
        {
            "request": request,
            "error": error,
            "qc_result": qc_result,
            "forward": forward,
            "reverse": reverse,
            "probe": probe or "",
            "template_sequence": template_sequence or "",
            "reference": reference,
        },
    )

