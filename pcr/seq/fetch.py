
import pysam

from pcr.components.variant import Variant
from pcr.seq.variant import apply_variant

from dataclasses import dataclass
from typing import Optional, Tuple, Dict, Any
import pysam
from Bio.Seq import Seq

@dataclass(frozen=True)
class GenomicRegion:
	chrom: str
	start: int  # 1-based inclusive
	end: int	# 1-based inclusive
	name: str

# --- 2. 유틸리티 함수: 특이도 강화를 위한 Mismatch 삽입 ---
def _apply_as_pcr_logic(
	target_base: str, 
	intensity: str = "strong"
) -> str:
	"""
	AS-PCR 프라이머의 특이도를 높이기 위해 3' 말단에서 n번째 위치에 Mismatch를 삽입합니다.
	"""
	# 3. 강한 Mismatch(Strong) 유도 (Purine <-> Purine, Pyrimidine <-> Pyrimidine 교체 등)
	mismatch_map = { 
		"strong": {'A': 'G', 'T': 'C', 'G': 'A', 'C': 'T'}, 
		"medium": {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}, 
		"weak": {'A': 'C', 'T': 'G', 'G': 'T', 'C': 'A'}, 
	}
	#mismatch_map = {
	#	"A": "G" if intensity == "strong" else "C",
	#	"G": "A" if intensity == "strong" else "T",
	#	"C": "T" if intensity == "strong" else "A",
	#	"T": "C" if intensity == "strong" else "G"
	#}
	# 해당 위치의 염기를 강제로 mismatch 염기로 교환
	return mismatch_map[intensity].get(target_base.upper(), "N")

def build_templates(
		fasta: pysam.FastaFile,
		region: GenomicRegion,
		ref_allele: str,
		alt_allele: str,
		max_amplicon_length: int,
		*,
		mismatch_pos: int = 3,
		max_primer_length: int = 30,
		ensure_pair_room: bool = True,
		debug: bool = True,
	) -> Dict[str, Any]:
	"""
	- genome reference window를 fetch한 뒤(reference_sequence)
	- target SNP(ref/alt) 반영 template 생성
	- mismatch_pos를 이용해 forward/reverse 모드별 mismatch index를 잡아
	  각 모드별 (ref/alt + mismatch) 템플릿까지 생성

	forward mismatch index = target_index - mismatch_pos
	reverse mismatch index = target_index + mismatch_pos
	"""

	(
		reference_sequence,
		t_start, t_end,
		target_index, target_end_index,
	) = fetch_template_sequence(
		fasta=fasta,
		region=region,
		max_amplicon_length=max_amplicon_length,
		max_primer_length=max_primer_length,
		ensure_pair_room=ensure_pair_room,
	)

	ref_allele = (ref_allele or "").upper().strip()
	alt_allele = (alt_allele or "").upper().strip()
	if len(ref_allele) != 1 or len(alt_allele) != 1:
		raise ValueError("ref_allele/alt_allele must be 1bp each.")

	# -------------------------
	# 1) target SNP 반영 템플릿
	# -------------------------
	var_ref = Variant(index=target_index, ref=ref_allele, alt=ref_allele)
	var_alt = Variant(index=target_index, ref=ref_allele, alt=alt_allele)

	ref_template = apply_variant(reference_sequence, var_ref, allele="alt")  # ref allele 유지
	alt_template = apply_variant(reference_sequence, var_alt, allele="alt")  # alt allele 반영

	# -------------------------
	# 2) mismatch index 계산 (요청대로 ± mismatch_pos)
	# -------------------------
	fw_mm_index = int(target_index) - int(mismatch_pos) + 1
	rv_mm_index = int(target_index) + int(mismatch_pos) - 1

	# mismatch index range 체크
	L = len(reference_sequence)
	def _check_mm_idx(idx: int, label: str) -> None:
		if idx < 0 or idx >= L:
			raise ValueError(
				f"{label} mismatch index out of range: idx={idx}, "
				f"target_index={target_index}, mismatch_pos={mismatch_pos}, len={L}"
			)

	_check_mm_idx(fw_mm_index, "forward")
	_check_mm_idx(rv_mm_index, "reverse")

	# -------------------------
	# 3) mismatch variant 생성 (reference base를 다른 염기로 강제 치환)
	#    - ref/alt 템플릿 각각에 동일 mismatch 적용 (시각화/QC 목적)
	# -------------------------
	# forward mismatch
	fw_ref_base = reference_sequence[fw_mm_index].upper()
	fw_mm_alt = _apply_as_pcr_logic(fw_ref_base)
	fw_mm_var = Variant(index=fw_mm_index, ref=fw_ref_base, alt=fw_mm_alt)

	# reverse mismatch
	rv_ref_base = reference_sequence[rv_mm_index].upper()
	rv_mm_alt = _apply_as_pcr_logic(rv_ref_base)
	rv_mm_var = Variant(index=rv_mm_index, ref=rv_ref_base, alt=rv_mm_alt)

	# apply mismatch onto ref/alt templates
	fw_ref_mm_template = apply_variant(ref_template, fw_mm_var, allele="alt")
	fw_alt_mm_template = apply_variant(alt_template, fw_mm_var, allele="alt")

	rv_ref_mm_template = apply_variant(ref_template, rv_mm_var, allele="alt")
	rv_alt_mm_template = apply_variant(alt_template, rv_mm_var, allele="alt")
	
	print(target_index)
	print(rv_mm_index)
	print(fw_mm_var)
	print(rv_mm_var)
	print(fw_ref_mm_template[target_index-2:])
	print(fw_alt_mm_template[target_index-2:])
	print(rv_ref_mm_template[target_index-2:])
	print(rv_alt_mm_template[target_index-2:])


	if debug:
		print(f"[DEBUG] Window: {region.chrom}:{t_start}-{t_end}")
		print(f"[DEBUG] Target Index in Template: {target_index}")
		print(f"[DEBUG] Base at Target (ref/alt): {ref_template[target_index]}/{alt_template[target_index]}")
		print(f"[DEBUG] Forward mismatch idx: {fw_mm_index} ref->{fw_ref_base} alt->{fw_mm_alt}")
		print(f"[DEBUG] Reverse mismatch idx: {rv_mm_index} ref->{rv_ref_base} alt->{rv_mm_alt}")

	return {
		# raw
		"reference_template_sequence": reference_sequence,

		# target index info (template 좌표)
		"target_start_index": target_index,
		"target_end_index": target_end_index,
		"template_start": t_start,
		"template_end": t_end,

		# base templates (SNP만 반영)
		"ref_template_sequence": ref_template,
		"alt_template_sequence": alt_template,

		# mode별 templates (SNP + mismatch까지 반영)
		"forward": {
			"mismatch_index": fw_mm_index,
			"mismatch_variant": fw_mm_var,
			"ref_template_sequence": ref_template,
			"alt_template_sequence": alt_template,
			"ref_mismatch_template_sequence": fw_ref_mm_template,
			"alt_mismatch_template_sequence": fw_alt_mm_template,
		},
		"reverse": {
			"mismatch_index": rv_mm_index,
			"mismatch_variant": rv_mm_var,
			"ref_template_sequence": ref_template,
			"alt_template_sequence": alt_template,
			"ref_mismatch_template_sequence": rv_ref_mm_template,
			"alt_mismatch_template_sequence": rv_alt_mm_template,
		},
	}
# --- 핵심 로직 함수 ---

def _get_contig_length(fasta: pysam.FastaFile, chrom: str) -> Optional[int]:
	try:
		return int(fasta.get_reference_length(chrom))
	except Exception:
		return None

def _clamp_window(start: int, end: int, contig_length: Optional[int]) -> Tuple[int, int]:
	start = max(1, start)
	if contig_length is not None:
		end = min(contig_length, end)
	if end < start:
		end = start
	return start, end

def compute_template_window_as_pcr(
	target_start: int,
	target_end: int,
	*,
	max_amplicon_length: Optional[int],
	max_primer_length: int = 30,
	contig_length: Optional[int] = None,
	ensure_pair_room: bool = True,
) -> Tuple[int, int]:
	"""
	AS-PCR 용 template window 계산 (1-based inclusive)
	Forward 고정을 위한 Upstream 공간과 Reverse 고정을 위한 Downstream 공간을 모두 확보합니다.
	"""
	flank_for_primer = max_primer_length - 1 

	# 1. Forward-fix 시: 타겟 왼쪽(flank) + 타겟 오른쪽(max_amplicon) 필요
	f_up = flank_for_primer
	f_down = (max_amplicon_length + max_primer_length) if (ensure_pair_room and max_amplicon_length) else flank_for_primer
	f_start, f_end = _clamp_window(target_start - f_up, target_end + f_down, contig_length)

	# 2. Reverse-fix 시: 타겟 오른쪽(flank) + 타겟 왼쪽(max_amplicon) 필요
	r_down = flank_for_primer
	r_up = (max_amplicon_length + max_primer_length) if (ensure_pair_room and max_amplicon_length) else flank_for_primer
	r_start, r_end = _clamp_window(target_start - r_up, target_end + r_down, contig_length)

	# 3. 두 경우를 모두 커버하는 UNION 영역 반환
	return min(f_start, r_start), max(f_end, r_end)

def fetch_template_sequence(
	fasta: pysam.FastaFile,
	region: GenomicRegion,
	*,
	max_amplicon_length: Optional[int],
	max_primer_length: int = 30,
	ensure_pair_room: bool = True,
) -> Tuple[str, int, int, int, int]:
	"""
	pysam을 이용해 서열을 가져오고 Primer3용 0-based 인덱스를 계산합니다.
	"""
	contig_len = _get_contig_length(fasta, region.chrom)

	template_start, template_end = compute_template_window_as_pcr(
		region.start, region.end,
		max_amplicon_length=max_amplicon_length,
		max_primer_length=max_primer_length,
		contig_length=contig_len,
		ensure_pair_room=ensure_pair_room,
	)

	# pysam.fetch(0-based start inclusive, 0-based end exclusive)
	template_sequence = fasta.fetch(
		region.chrom, template_start - 1, template_end
	).upper()

	# ✅ [중요] 타겟 인덱스 계산 (0-based)
	# 예: region.start가 10이고 template_start가 1이면, 인덱스는 9입니다.
	target_start_index = region.start - template_start
	target_len = region.end - region.start + 1
	target_end_index = target_start_index + target_len - 1

	return template_sequence, template_start, template_end, target_start_index, target_end_index

def build_ref_alt_templates(
	fasta: pysam.FastaFile,
	region: GenomicRegion,
	ref_allele: str,
	alt_allele: str,
	max_amplicon_length: int,
	max_primer_length: int = 30,
	ensure_pair_room: bool = True,
	debug: bool = True,
) -> Dict[str, object]:
		
	(
		reference_sequence,
		t_start, t_end,
		t_start_idx, t_end_idx
	) = fetch_template_sequence(
		fasta=fasta, region=region,
		max_amplicon_length=max_amplicon_length,
		max_primer_length=max_primer_length,
		ensure_pair_room=ensure_pair_room
	)

	# Variant 적용
	var_ref = Variant(index=t_start_idx, ref=ref_allele, alt=ref_allele)
	var_alt = Variant(index=t_start_idx, ref=ref_allele, alt=alt_allele)

	ref_template_sequence = apply_variant(reference_sequence, var_ref, allele="alt")
	alt_template_sequence = apply_variant(reference_sequence, var_alt, allele="alt")

	if debug:
		print(f"[DEBUG] Window: {region.chrom}:{t_start}-{t_end}")
		print(f"[DEBUG] Target Index in Template: {t_start_idx}")
		print(f"[DEBUG] Base at Index (Ref/Alt): {ref_template_sequence[t_start_idx]}/{alt_template_sequence[t_start_idx]}")

	return {
		"reference_template_sequence": reference_sequence,
		"ref_template_sequence": ref_template_sequence,
		"alt_template_sequence": alt_template_sequence,
		"target_start_index": t_start_idx,
		"target_end_index": t_end_idx,
		"template_start": t_start,
		"template_end": t_end,
	}