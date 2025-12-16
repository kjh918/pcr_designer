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


# @router.post("/excel")
# async def qc_excel_run(
#     reference: str = Form("hg38"),
#     file: UploadFile = File(None),
# ):
#     """
#     Excel/TSV 파일로 여러 primer pair를 일괄 QC.
#     기대 컬럼:
#       - Forward_Primer (이름, optional)
#       - Forward_seq   (필수)
#       - Reverse_Primer
#       - Reverse_seq
#       - Probe_seq     (optional)
#     """
#     # 파일 읽기
#     content = await file.read()
#     buf = io.BytesIO(content)

#     # 확장자에 따라 파싱
#     filename = file.filename or ""
#     filename_lower = filename.lower()

#     if filename_lower.endswith((".xlsx", ".xls")):
#         df = pd.read_excel(buf)
#     else:
#         # 기본은 TSV로 가정 (원래 스크립트도 tsv 썼으니까)
#         df = pd.read_csv(buf, sep="\t")

#     # BLAST DB
#     ref_cfg = settings.references[reference]
#     blast_db = str(ref_cfg.blast)

#     qc_thresholds = QCThresholds()

#     result_rows: list[dict[str, Any]] = []

#     for _, row in df.iterrows():
#         f_name = (
#             row.get("Forward_Primer")
#             if "Forward_Primer" in df.columns
#             else "FORWARD"
#         )
#         r_name = (
#             row.get("Reverse_Primer")
#             if "Reverse_Primer" in df.columns
#             else "REVERSE"
#         )

#         f_seq = str(row["Forward_seq"]).strip()
#         r_seq = str(row["Reverse_seq"]).strip()
#         p_seq = (
#             str(row["Probe_seq"]).strip()
#             if "Probe_seq" in df.columns and pd.notna(row["Probe_seq"])
#             else None
#         )

#         # 1) Thermo QC
#         amp = make_amplicon_for_qc(
#             forward_seq=f_seq,
#             reverse_seq=r_seq,
#             probe_seq=p_seq,
#             template_seq=None,
#         )
#         thermo_total, _ = evaluate_amplicons(
#             genomic_id=str(f_name),
#             amplicons=[amp],
#             qc_thresholds=qc_thresholds,
#         )
#         thermo_row = thermo_total[0]

#         # 2) BLAST QC
#         blast_qc = blast_qc_for_primer_pair(
#             f_name=str(f_name),
#             f_seq=f_seq,
#             r_name=str(r_name),
#             r_seq=r_seq,
#             db=blast_db,
#             probe_name="PROBE" if p_seq else None,
#             probe_seq=p_seq,
#         )

#         # 3) 원본 input + QC 결과 합치기
#         out_row = {
#             # 원본 input 보존
#             "Forward_Primer": f_name,
#             "Forward_seq": f_seq,
#             "Reverse_Primer": r_name,
#             "Reverse_seq": r_seq,
#             "Probe_seq": p_seq or "",
#             # Thermo (evaluate_amplicons 결과)
#             **thermo_row,
#             # BLAST
#             "BLAST_F_HITS": blast_qc["f_hits"],
#             "BLAST_R_HITS": blast_qc["r_hits"],
#             "BLAST_NEARBY_AMP_COUNT": blast_qc["nearby_count"],
#             "BLAST_MIN_AMP_SIZE": blast_qc["min_amplicon_size"],
#             "BLAST_AMP_DETAIL": ";".join(blast_qc["amplicon_details"])
#             if blast_qc["amplicon_details"]
#             else "",
#             "QC_BLAST_HIT": blast_qc["qc_blast_hit"],
#             "QC_BLAST_AMP": blast_qc["qc_blast_amplicon"],
#             "QC_PROBE_IN_AMP": blast_qc.get("qc_probe_in_amplicon", "-"),
#         }

#         result_rows.append(out_row)

#     result_df = pd.DataFrame(result_rows)

#     # 엑셀로 만들어서 반환
#     out_buf = io.BytesIO()
#     with pd.ExcelWriter(out_buf, engine="xlsxwriter") as writer:
#         result_df.to_excel(writer, index=False, sheet_name="QC_Result")

#     out_buf.seek(0)

#     return StreamingResponse(
#         out_buf,
#         media_type=(
#             "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
#         ),
#         headers={
#             "Content-Disposition": 'attachment; filename="qc_result.xlsx"'
#         },
#     )
