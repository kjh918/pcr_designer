# app/routers/design_aspcr.py
from __future__ import annotations

from typing import Any, Dict, List, Optional
import traceback

from fastapi import APIRouter, Request, Depends, Form, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from datetime import datetime

from app.routers.design_common import (
	CommonDesignForm,
	build_common_kwargs,
	parse_regions_from_form,
	init_context,
)


from pcr.config.runtime import (
	get_fasta_handle,
	build_pcr_params_from_web,
	build_qc_params_from_web,
	resolve_pcr_params,
	merge_dict,
)


from pcr.seq.fetch import GenomicRegion, build_ref_alt_templates
from pcr.designers.as_pcr import AsPcrDesigner
from pcr.pipelines.as_pcr import run_as_pcr_pipeline


router = APIRouter(prefix="/design", tags=["design"])
templates = Jinja2Templates(directory="app/templates")


@router.post("/as-pcr", response_class=HTMLResponse, name="design_aspcr_from_form")
async def design_aspcr_from_form(
	request: Request,
	f: CommonDesignForm = Depends(CommonDesignForm.as_form),

	# ✅ form_aspcr.html : name="ref", name="alt"
	ref_allele: str = Form("", alias="ref"),
	alt_allele: str = Form("", alias="alt"),

	# (선택) 후보 개수 제한
	n_primers: Optional[int] = Form(None),

	# (선택) forward 길이 범위/ reverse 후보 수 (폼에서 받으면 연결)
	forward_min_len: Optional[int] = Form(None),
	forward_max_len: Optional[int] = Form(None),
	n_reverse: Optional[int] = Form(None),
):
	context = init_context(request, f, assay="as-pcr")

	try:
		ref_allele = (ref_allele or "").strip()
		alt_allele = (alt_allele or "").strip()

		if not ref_allele or not alt_allele:
			raise HTTPException(status_code=400, detail="AS-PCR은 ref/alt 입력이 필요합니다.")
		if len(ref_allele) != 1 or len(alt_allele) != 1:
			raise HTTPException(status_code=400, detail="AS-PCR(ref/alt)은 현재 1bp(SNP)만 지원합니다.")

		# FASTA handle
		fasta = get_fasta_handle(f.reference)

		# region(s)
		regions = await parse_regions_from_form(f)
		if not regions:
			raise HTTPException(status_code=400, detail="입력 region이 없습니다.")
		print(regions)
		# 공통 kwargs
		common_kwargs = build_common_kwargs(f)

		# -----------------------------
		# ✅ 1) request-scope PCR/QC config 생성 (qpcr 방식)
		# -----------------------------
		pcr_overrides: Dict[str, Any] = {
			"primer_kwargs": {
				"min_amplicon_length": common_kwargs.get("min_amplicon_length"),
				"max_amplicon_length": common_kwargs.get("max_amplicon_length"),
				"n_primers": n_primers if n_primers is not None else common_kwargs.get("n_primers"),
				"primer3_global_args": {
					# AS-PCR primer3 인자(필요하면 여기에 매핑)
					# 예: "PRIMER_OPT_SIZE": ...
				},
			},
			"bisulfite": {"run": False},
		}
		pcr_cfg = build_pcr_params_from_web(pcr_overrides)

		qc_overrides: Dict[str, Any] = {
			"MIN_AMP_BP": common_kwargs.get("min_amplicon_length"),
			"MAX_AMP_BP": common_kwargs.get("max_amplicon_length"),
			# thermo dG, dimer 등 폼 연결되면 여기에 추가
		}
		qc_params = build_qc_params_from_web(qc_overrides)

		resolved = resolve_pcr_params(
			pcr_cfg=pcr_cfg,
			min_amplicon_length=common_kwargs.get("min_amplicon_length"),
			max_amplicon_length=common_kwargs.get("max_amplicon_length"),
			n_primers=n_primers if n_primers is not None else common_kwargs.get("n_primers"),
			bisulfite=False,
		)

		primer3_global_args = merge_dict(
			base=pcr_cfg.primer_kwargs.primer3_global_args,
			override=None,
		)

		if f.mode == "single":
			r = regions[0]
			gr = GenomicRegion(chrom=r.chrom, start=r.start, end=r.end, name=f.name)

			result = run_as_pcr_pipeline(
				genomic_id=f.name,
				region=gr,
				fasta=fasta,
				ref_allele=ref_allele,
				alt_allele=alt_allele,
				pcr_cfg=pcr_cfg,
				qc_params=qc_params,
				min_amplicon_length=resolved.min_amplicon_length,
				max_amplicon_length=resolved.max_amplicon_length,
				forward_max_len=30,
				forward_min_len=15,
				mismatch_pos=3,
				n_primers=resolved.n_primers,
				n_reverse=100,
				primer3_global_args=primer3_global_args,
			)

			total_df = result.total_df
			filtered_df = result.filtered_df
			
			context["single_result"] = {
				"region": r,
				"total_count": len(total_df),
				"filtered_count": len(filtered_df),
			}
			context["single_total_amplicons"] = total_df.to_dict(orient="records")
			context["single_filtered_amplicons"] = filtered_df.to_dict(orient="records")
			context["export_meta"] = {
				# ---- 기본 정보 ----
				"assay": "as-pcr",
				"timestamp": datetime.now().isoformat(timespec="seconds"),

				# ---- 입력 정보 ----
				"reference": f.reference,
				"region": {
					"chrom": r.chrom,
					"start": r.start,
					"end": r.end,
					"name": r.name,
				},
				"ref_allele": ref_allele,
				"alt_allele": alt_allele,

				# ---- 결과 요약 ----
				"total_count": len(total_df),
				"filtered_count": len(filtered_df),

				# ---- resolved PCR params (재현성 핵심) ----
				"pcr_params": {
					"min_amplicon_length": resolved.min_amplicon_length,
					"max_amplicon_length": resolved.max_amplicon_length,
					"n_primers": resolved.n_primers,
					"forward_min_len": forward_min_len,
					"forward_max_len": forward_max_len,
					"n_reverse": n_reverse,
				}
			}



		else:
			multi_results: List[Dict[str, Any]] = []

			for r in regions:
				gr = GenomicRegion(chrom=r.chrom, start=r.start, end=r.end, name=r.name)

				result = run_as_pcr_pipeline(
					genomic_id=f.name,
					region=gr,
					fasta=fasta,
					ref_allele=ref_allele,
					alt_allele=alt_allele,
					pcr_cfg=pcr_cfg,
					qc_params=qc_params,
					min_amplicon_length=resolved.min_amplicon_length,
					max_amplicon_length=resolved.max_amplicon_length,
					n_primers=resolved.n_primers,
					n_reverse=resolved.n_primers,
					forward_min_len=18,
					forward_max_len=28,
					primer3_global_args=primer3_global_args,
				)

				total_df = result.total_df
				filtered_df = result.filtered_df

				multi_results.append(
					dict(
						region=r,
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
