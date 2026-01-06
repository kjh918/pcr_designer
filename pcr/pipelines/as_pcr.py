# pcr/pipelines/as_pcr.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, List

import pandas as pd

from pcr.seq.fetch import GenomicRegion, fetch_template_sequence, build_templates
from pcr.config.schema.qc import QCParams
from pcr.designers.as_pcr import AsPcrDesigner
from pcr.pipelines.base import run_pipeline_from_amplicons, PipelineResult
import primer3

def _pick_template_for_run(
	ref_alt_templates_dict: Dict[str, Any],
	*,
	template_type: str,
	fixed_primer: str,  # "forward" or "reverse"
) -> str:
	"""
	template_type(wt/alt/wt_mm/alt_mm) + fixed_primer(forward/reverse)에 따라
	AsPcrDesigner에 넣을 template_sequence를 선택한다.
	"""
	template_type = (template_type or "").lower()
	fixed_primer = (fixed_primer or "").lower()

	if template_type == "wt":
		return ref_alt_templates_dict["ref_template_sequence"]
	if template_type == "alt":
		return ref_alt_templates_dict["alt_template_sequence"]
	if template_type == "wt_mm":
		return ref_alt_templates_dict[fixed_primer]["ref_mismatch_template_sequence"]
	if template_type == "alt_mm":
		return ref_alt_templates_dict[fixed_primer]["alt_mismatch_template_sequence"]

	raise ValueError(f"Unknown template_type: {template_type}")

def run_as_pcr_pipeline(
	*,
	genomic_id: str,
	region: GenomicRegion,
	fasta: pysam.FastaFile,
	ref_allele: str,
	alt_allele: str,
	mismatch_pos:int,
	pcr_cfg: Any,  # (지금은 안 쓰지만 유지 가능)
	qc_params: QCParams,  # ✅ qc_cfg → qc_params 로 명확히
	min_amplicon_length: Optional[int],
	max_amplicon_length: Optional[int],
	n_primers: Optional[int],
	n_reverse: Optional[int],
	forward_max_len: Optional[int],
	forward_min_len: Optional[int],
	primer3_global_args: Optional[Dict],
) -> AsPipelineResult:
	"""
	AS-PCR 후보 생성 + Thermo QC + BLAST QC까지 수행 후
	total_df / filtered_df 반환.

	filtered_df 기준:
	  - WT/ALT 모두 THERMO_PASS == O
	  - WT/ALT 모두 BLAST_PASS == O
	"""

	ref_alt_templates_dict = build_templates(
		region=region,
		fasta=fasta,
		ref_allele=ref_allele,
		alt_allele=alt_allele,
		max_amplicon_length=max_amplicon_length,
		mismatch_pos=mismatch_pos,
	)

	reference_window = ref_alt_templates_dict["reference_template_sequence"]  # "진짜 ref"
	ref_tpl = ref_alt_templates_dict["ref_template_sequence"]				# SNP(ref)
	alt_tpl = ref_alt_templates_dict["alt_template_sequence"]				# SNP(alt)

	# ---- mode templates (SNP + mismatch) ----
	fw_ref_mm_tpl = ref_alt_templates_dict["forward"]["ref_mismatch_template_sequence"]
	fw_alt_mm_tpl = ref_alt_templates_dict["forward"]["alt_mismatch_template_sequence"]

	rv_ref_mm_tpl = ref_alt_templates_dict["reverse"]["ref_mismatch_template_sequence"]
	rv_alt_mm_tpl = ref_alt_templates_dict["reverse"]["alt_mismatch_template_sequence"]

	# ---- indices ----
	t_start_idx = int(ref_alt_templates_dict["target_start_index"])
	t_end_idx = int(ref_alt_templates_dict["target_end_index"])

	chrom = str(getattr(ref_alt_templates_dict.get("region", None), "chrom", ref_alt_templates_dict.get("chrom", "")) or "")
	start = int(ref_alt_templates_dict.get("template_start", 0))
	end = int(ref_alt_templates_dict.get("template_end", len(reference_window) - 1))
	ref_allele = (ref_alt_templates_dict.get("ref_allele") or "").upper()
	alt_allele = (ref_alt_templates_dict.get("alt_allele") or "").upper()
	
	total_amplicons: List[Amplicon] = []

	# ======================================================================
	# 1) forward fixed: primer3는 WT(ref template) 기준으로 1번만 실행
	#	-> 좌표로 wt/alt/wt_mm/alt_mm 4종을 "세트"로 만든다
	# ======================================================================
	fw_designer = AsPcrDesigner(
		reference_template_sequence=ref_tpl,	   # wt
		alt_template_sequence=alt_tpl,			 # alt
		ref_mm_template_sequence=fw_ref_mm_tpl,	# wt_mm (forward mismatch)
		alt_mm_template_sequence=fw_alt_mm_tpl,	# alt_mm (forward mismatch)

		target_start_index=t_start_idx,
		target_end_index=t_end_idx,
		target_index=t_start_idx,

		ref_allele=ref_allele if ref_allele else ref_tpl[t_start_idx],
		alt_allele=alt_allele if alt_allele else alt_tpl[t_start_idx],

		chrom=chrom,
		start=start,
		end=end,

		mismatch_pos=mismatch_pos,
		fixed_prime="forward",

		min_amplicon_length=min_amplicon_length,
		max_amplicon_length=max_amplicon_length,
		n_primers=int(n_reverse) if n_reverse is not None else 100,
		primer3_global_args=primer3_global_args,
	)

	fw_sets = fw_designer.design_sets()   # ✅ set 단위 생성
	total_amplicons.extend(fw_designer.amplicon_list)  # ✅ flat

	# ======================================================================
	# 2) reverse fixed: primer3 1번 실행 -> 4종 세트 생성
	# ======================================================================
	rv_designer = AsPcrDesigner(
		reference_template_sequence=ref_tpl,
		alt_template_sequence=alt_tpl,
		ref_mm_template_sequence=rv_ref_mm_tpl,	# wt_mm (reverse mismatch)
		alt_mm_template_sequence=rv_alt_mm_tpl,	# alt_mm (reverse mismatch)

		target_start_index=t_start_idx,
		target_end_index=t_end_idx,
		target_index=t_start_idx,

		ref_allele=ref_allele if ref_allele else ref_tpl[t_start_idx],
		alt_allele=alt_allele if alt_allele else alt_tpl[t_start_idx],

		chrom=chrom,
		start=start,
		end=end,

		mismatch_pos=mismatch_pos,
		fixed_prime="reverse",

		min_amplicon_length=min_amplicon_length,
		max_amplicon_length=max_amplicon_length,
		n_primers=int(n_reverse) if n_reverse is not None else 100,
		primer3_global_args=primer3_global_args,
	)

	rv_sets = rv_designer.design_sets()
	total_amplicons.extend(rv_designer.amplicon_list)

	return run_pipeline_from_amplicons(
		genomic_id=genomic_id,
		qc_params=qc_params,
		assay="aspcr",
		amplicon_list=total_amplicons)
	# ✅ QCThresholds 생성 제거 → QCParams 그대로 주입
	
