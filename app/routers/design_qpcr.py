# app/routers/design_qpcr.py
from __future__ import annotations

from typing import Any, Dict, List
import traceback
from datetime import datetime

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

from pcr.config.runtime import (
	get_fasta_handle,
	build_pcr_params_from_web,
	build_qc_params_from_web,
	resolve_pcr_params,
	merge_dict,
)

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

		# -----------------------------
		# ✅ 1) 웹 입력 기반 request-scope PCR/QC config 생성
		# -----------------------------
		pcr_overrides: Dict[str, Any] = {
			"primer_kwargs": {
				# 공통 폼 값 (있으면 덮어씀)
				"min_amplicon_length": common_kwargs.get("min_amplicon_length"),
				"max_amplicon_length": common_kwargs.get("max_amplicon_length"),
				"n_primers": common_kwargs.get("n_primers"),
				# primer3 args는 (추후 폼 연결되면 여기에 추가)
				# "primer3_global_args": {...}
			},
			"probe_kwargs": {
				"n_probes": effective_n_probes,
				"primer3_global_args": {
					# probe 관련 입력을 primer3 internal probe 키로 매핑
					# (프로젝트에서 쓰는 키 네이밍이 다르면 여기만 조정)
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
			# diff tm
			"PROBE_MIN_DIFF_TM": min_primer_probe_tm_diff,
			"PROBE_MAX_DIFF_TM": max_primer_probe_tm_diff,

			# amplicon QC filter도 동일 범위로 동기화(혼선 방지)
			"MIN_AMP_BP": common_kwargs.get("min_amplicon_length"),
			"MAX_AMP_BP": common_kwargs.get("max_amplicon_length"),
		}

		qc_cfg = build_qc_params_from_web(qc_overrides)

		# -----------------------------
		# ✅ 2) resolve는 merged pcr_cfg 기준으로!
		# -----------------------------
		resolved = resolve_pcr_params(
			pcr_cfg=pcr_cfg,
			min_amplicon_length=common_kwargs.get("min_amplicon_length"),
			max_amplicon_length=common_kwargs.get("max_amplicon_length"),
			n_probes=effective_n_probes,
			n_primers=common_kwargs.get("n_primers"),
			bisulfite=False,
		)

		# primer3 args는 config 기본 + override merge
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
			r = regions[0]
			gr = GenomicRegion(chrom=region.chrom, start=region.start, end=region.end, name=region.name)

			result = run_qpcr(
				region=gr,
				fasta=fasta,
				pcr_cfg=pcr_cfg,   # ✅ merged
				qc_params=qc_cfg,	 # ✅ merged
				min_amplicon_length=resolved.min_amplicon_length,
				max_amplicon_length=resolved.max_amplicon_length,
				n_probes=resolved.n_probes,
				n_primers=resolved.n_primers,
				primer3_global_args=primer3_global_args,
				probe_primer3_global_args=probe_primer3_global_args,
			)

			total_df = result.total_df
			filtered_df = result.filtered_df
			print(total_df.columns)
			#x
			context["single_result"] = {
				"region": r,
				"total_count": len(total_df),
				"filtered_count": len(filtered_df),
			}
			context["single_total_amplicons"] = total_df.to_dict(orient="records")
			context["single_filtered_amplicons"] = filtered_df.to_dict(orient="records")
			context["export_meta"] = {
				# ---- 기본 정보 ----
				"assay": "qpcr",
				"timestamp": datetime.now().isoformat(timespec="seconds"),

				# ---- 입력 정보 ----
				"reference": f.reference,
				"region": {
					"chrom": r.chrom,
					"start": r.start,
					"end": r.end,
					"name": r.name,
				},
				# ---- 결과 요약 ----
				"total_count": len(total_df),
				"filtered_count": len(filtered_df),

				# ---- resolved PCR params (재현성 핵심) ----
				"pcr_params": {
					"min_amplicon_length": resolved.min_amplicon_length,
					"max_amplicon_length": resolved.max_amplicon_length,
					"n_primers": resolved.n_primers,
				}
			}

		else:
			multi_results: List[Dict[str, Any]] = []
			for region in regions:
				gr = GenomicRegion(chrom=region.chrom, start=region.start, end=region.end, name=region.name)
				
				result = run_qpcr(
					region=gr,
					fasta=fasta,
					pcr_cfg=pcr_cfg,  # ✅ merged
					qc_params=qc_cfg,	# ✅ merged
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
