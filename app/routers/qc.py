# app/routers/qc.py
from __future__ import annotations

from fastapi import APIRouter, Request, Form
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from pcr.config.settings import settings
from pcr.config.runtime import build_qc_params_from_web

from pcr.qc.evaluate import run_thermo_qc, run_blast_qc
from pcr.qc.factories import make_amplicons_for_qc_from_blast

templates = Jinja2Templates(directory="app/templates")

router = APIRouter(prefix="/qc", tags=["qc"])


@router.get("/", response_class=HTMLResponse)
async def qc_page(request: Request, error: str | None = None):
	"""
	QC Only 화면 (초기 상태)
	"""
	return templates.TemplateResponse(
		"qc.html",
		{
			"request": request,
			"error": error,

			# ✅ CHANGED: qpcr_single과 동일한 키로 렌더되게
			"single_result": None,
			"single_total_amplicons": None,
			"single_filtered_amplicons": None,
			"export_meta": None,

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
	template_sequence: str | None = Form(None),  # 화면 호환용(현재 blast-template 사용 시 미사용)
	# --- reference ---
	reference: str = Form("hg38"),
	# --- QC override ---
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
	✅ QC Only도 qpcr_single과 동일한 output shape로 렌더:
	  - single_result
	  - single_total_amplicons
	  - single_filtered_amplicons
	  - export_meta
	"""
	error: str | None = None

	# ✅ CHANGED: qpcr_single 동일 키
	single_result = None
	single_total_amplicons: list[dict] | None = None
	single_filtered_amplicons: list[dict] | None = None
	export_meta: dict | None = None

	try:
		# -------------------------------------------------
		# 1) QCParams 생성 (settings + web override)
		# -------------------------------------------------
		qc_params = build_qc_params_from_web(
			{
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
			}
		)

		# -------------------------------------------------
		# 2) reference 기반 blast db / fasta
		# -------------------------------------------------
		if reference not in settings.references:
			raise ValueError(f"Unknown reference: {reference}")

		ref_cfg = settings.references[reference]
		blast_db = str(ref_cfg.blast)
		fasta_path = str(ref_cfg.fasta)

		# -------------------------------------------------
		# 3) BLAST → amplicon 후보들 생성
		# -------------------------------------------------
		blast_qc = run_blast_qc(
			f_name="FORWARD",
			f_seq=forward,
			r_name="REVERSE",
			r_seq=reverse,
			db=blast_db,
			fasta=fasta_path,
			qc_params=qc_params,
			probe_name="PROBE" if probe else None,
			probe_seq=probe,
		)

		if blast_qc.get("blast_error") is True:
			raise RuntimeError("BLAST error")

		# -------------------------------------------------
		# 4) BLAST PASS된 후보들만 Amplicon으로 변환 → Thermo QC 입력
		# -------------------------------------------------
		amplicons = make_amplicons_for_qc_from_blast(
			blast_result=blast_qc,
			forward_seq=forward,
			reverse_seq=reverse,
			probe_seq=probe,
		)

		if not amplicons:
			# ✅ CHANGED: qpcr_single shape는 유지하되, count=0으로 표시
			single_result = {
				"region": {
					"chrom": "QC_ONLY",
					"start": 0,
					"end": 0,
				},
				"total_count": 0,
				"filtered_count": 0,
			}
			single_total_amplicons = []
			single_filtered_amplicons = []
			export_meta = {
				"kind": "qc_only",
				"reference": reference,
				"blast_db": blast_db,
				"fasta": fasta_path,
				"message": "BLAST PASS된 amplicon 후보가 없습니다.",
			}

			return templates.TemplateResponse(
				"qc.html",
				{
					"request": request,
					"error": None,
					"single_result": single_result,
					"single_total_amplicons": single_total_amplicons,
					"single_filtered_amplicons": single_filtered_amplicons,
					"export_meta": export_meta,
					"forward": forward,
					"reverse": reverse,
					"probe": probe or "",
					"template_sequence": template_sequence or "",
					"reference": reference,
				},
			)

		# -------------------------------------------------
		# 5) Thermo QC
		# -------------------------------------------------
		total_rows, filtered_rows = run_thermo_qc(
			genomic_id="QC_ONLY",
			amplicons=amplicons,
			qc_params=qc_params,
		)

		# -------------------------------------------------
		# 6) qpcr_single과 동일 shape로 매핑
		# -------------------------------------------------
		single_total_amplicons = total_rows or []
		single_filtered_amplicons = filtered_rows or []

		# ✅ CHANGED: QC Only region은 blast 결과로 "대표 region" 하나 만들어주기
		# - 기존 results_single.html이 region.start/end를 찍으니까 최소한 숫자는 넣어야 함
		# - blast_qc.result.filtered_amplicons[0]의 reference 좌표를 대표로 사용
		rep_chrom = "QC_ONLY"
		rep_start = 0
		rep_end = 0

		blast_res = blast_qc.get("result") or {}
		filt_amps = blast_res.get("filtered_amplicons") or []
		print(single_total_amplicons)
		if isinstance(filt_amps, list) and filt_amps:
			a0 = filt_amps[0]
			rep_chrom = (a0.get("genomic") or {}).get("chrom") or rep_chrom
			ref0 = a0.get("reference") or {}
			rep_start = int(ref0.get("start_1b") or 0)
			rep_end = int(ref0.get("end_1b") or 0)

		single_result = {
			"region": {
				"chrom": rep_chrom,
				"start": rep_start,
				"end": rep_end,
			},
			"total_count": len(single_total_amplicons),
			"filtered_count": len(single_filtered_amplicons),
		}

		# ✅ CHANGED: export_meta에 BLAST 후보 요약(좌표/말단 mismatch 등) 넣기
		amp_summaries = []
		for i, a in enumerate(filt_amps[:200], start=1):  # 너무 커지는 것 방지
			ref = a.get("reference") or {}
			amp = a.get("amplicon") or {}
			bind = a.get("binding") or {}
			f_bind = (bind.get("forward") or {})
			r_bind = (bind.get("reverse") or {})
			amp_summaries.append(
				{
					"idx": i,
					"chrom": (a.get("genomic") or {}).get("chrom"),
					"ref_start_1b": ref.get("start_1b"),
					"ref_end_1b": ref.get("end_1b"),
					"amp_start_1b": amp.get("start_1b"),
					"amp_end_1b": amp.get("end_1b"),
					"amp_len": amp.get("length_bp"),
					"forward_unmatched_3p": len(f_bind.get("unmatched_3p_indices") or []),
					"reverse_unmatched_3p": len(r_bind.get("unmatched_3p_indices") or []),
					"probe_in_amplicon": bind.get("probe_in_amplicon"),
				}
			)

		export_meta = {
			"kind": "qc_only",
			"reference": reference,
			"blast_db": blast_db,
			"fasta": fasta_path,
			"input": {
				"forward": forward,
				"reverse": reverse,
				"probe": probe or "",
			},
			"blast_filtered_amplicon_count": len(filt_amps),
			"blast_amplicon_summaries": amp_summaries,
		}

	except Exception as e:
		error = f"QC 수행 중 오류: {e}"

	return templates.TemplateResponse(
		"qc.html",
		{
			"request": request,
			"error": error,

			# ✅ CHANGED: qpcr_single과 동일 키
			"single_result": single_result,
			"single_total_amplicons": single_total_amplicons,
			"single_filtered_amplicons": single_filtered_amplicons,
			"export_meta": export_meta,

			# form 값 유지
			"forward": forward,
			"reverse": reverse,
			"probe": probe or "",
			"template_sequence": template_sequence or "",
			"reference": reference,
		},
	)
