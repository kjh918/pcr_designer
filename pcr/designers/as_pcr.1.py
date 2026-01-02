from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple
import primer3
from Bio.Seq import Seq

from pcr.designers.base import BasePrimerDesigner
from pcr.components import Amplicon
from pcr.components.primer import Primer

# ----------------------------------------------------------------------
# Helpers (Modified & Integrated)
# ----------------------------------------------------------------------

def _replace_3prime_base(primer_5to3: str, base: str, strand: str = "forward") -> str:
	"""
	Replaces the 3' end of the primer with the allele base.
	For Reverse primers, uses the complement of the allele.
	"""
	s = list(primer_5to3.upper())
	if not s: return ""
		
	if strand == "forward":
		s[-1] = base.upper()
	else:
		# Reverse primer binds to the sense strand, so it needs the complement
		s[-1] = str(Seq(base).complement()).upper()
	return "".join(s)

def _apply_as_pcr_logic(
	primer_5to3: str,
	target_base: str,
	strand: str,
	mismatch_pos: int = 3,
	intensity: str = "strong",
) -> str:
	"""
	Applies AS-PCR logic: 3' end replacement + internal mismatch for specificity.
	"""
	mismatch_map = {
		"strong": {'A': 'G', 'T': 'C', 'G': 'A', 'C': 'T'},
		"medium": {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'},
		"weak":   {'A': 'C', 'T': 'G', 'G': 'T', 'C': 'A'},
	}

	# 1. Set the 3' end first
	primer_seq = _replace_3prime_base(primer_5to3, target_base, strand)
	s = list(primer_seq)
		
	# 2. Apply internal mismatch
	idx = -int(mismatch_pos)
	if abs(idx) <= len(s):
		original = s[idx]
		new_base = mismatch_map.get(intensity.lower(), mismatch_map["strong"]).get(original, "A")
		if new_base == original: # Safety check
			for cand in ("A", "C", "G", "T"):
				if cand != original:
					new_base = cand
					break
		s[idx] = new_base
		
	return "".join(s)

def _patch_template_for_pair(
		template_5to3: str,
		*,
		left_seq_5to3: str,
		left_start: int,
		left_end: int,
		right_seq_5to3: str,
		right_start: int,
		right_end: int,
	) -> str:
	"""Overwrites the template with actual primer sequences (Right is RC'd)."""
	t = list(template_5to3.upper())
	left_seq = left_seq_5to3.upper()
	right_on_template = str(Seq(right_seq_5to3.upper()).reverse_complement())

	t[left_start : left_end + 1] = list(left_seq)
	t[right_start : right_end + 1] = list(right_on_template)
	return "".join(t)

def _mk_primer(
		*,
		template_sequence: str,
		reference_template_sequence: str,
		sequence: str,
		strand: str,
		primer_type: str,
		target_start_index: int,
		target_end_index: int,
		binding_start_index: int,
		binding_end_index: int,
	) -> Primer:
	return Primer(
		template_sequence=template_sequence,
		reference_template_sequence=reference_template_sequence,
		sequence=sequence,
		strand=strand,
		primer_type=primer_type,
		target_start_index=target_start_index,
		target_end_index=target_end_index,
		binding_start_index=binding_start_index,
		binding_end_index=binding_end_index,
	)

def _pos_to_binding_strict(pos: Any) -> Tuple[int, int]:
	if not pos or len(pos) != 2:
		raise ValueError(f"Invalid position: {pos}")
	start, length = int(pos[0]), int(pos[1])
	return start, start + length - 1

# ----------------------------------------------------------------------
# AsPcrDesigner
# ----------------------------------------------------------------------

class AsPcrDesigner(BasePrimerDesigner):
	def __init__(
		self,
		template_sequence: str,
		reference_template_sequence: str,
		target_start_index: int,
		target_end_index: int,
		target_index: int,
		ref_allele: str,
		alt_allele: str,
		**kwargs
	) -> None:
		self.target_index = int(target_index)
		self.ref_allele = ref_allele.upper()
		self.alt_allele = alt_allele.upper()
		super().__init__(
			template_sequence=template_sequence,
			reference_template_sequence=reference_template_sequence,
			target_start_index=target_start_index,
			target_end_index=target_end_index,
			**kwargs
		)
		
	def reset(self) -> None:
		"""Primer3 설정 인자를 초기화합니다."""
		self.primer3_seq_args = {
			'SEQUENCE_TEMPLATE': self.template_sequence,
		}
		# 기본 global 인자 설정 (필요에 따라 수정 가능)
		self.primer3_global_args = {
			'PRIMER_OPT_SIZE': 20,
			'PRIMER_MIN_SIZE': 18,
			'PRIMER_MAX_SIZE': 25,
			'PRIMER_OPT_TM': 55.0,
			'PRIMER_MIN_TM': 45.0,
			'PRIMER_MAX_TM': 65.0,
			'PRIMER_OPT_GC': 45.0,
			'PRIMER_MIN_GC': 20.0,
			'PRIMER_MAX_GC': 85.0,
			'PRIMER_MAX_POLY_X': 5,
			'PRIMER_SALT_MONOVALENT': 50.0,
			'PRIMER_DNA_CONC': 50.0,
		}

	def _configure_primer_forward_fix(self) -> None:
		"""
		Forward 프라이머의 3' end를 target_index에 고정하여 
		해당 지점에서 증폭이 시작되도록 설정합니다.
		"""
		super()._configure_primer_common()

		# 1. AS-PCR은 변이 지점 자체를 프라이머가 물고 있어야 하므로 
		# 기존에 설정된 SEQUENCE_TARGET(회피 영역 등)이 있다면 제거합니다.
		self.primer3_seq_args.pop("SEQUENCE_TARGET", None)

		# 2. ✅ LEFT(Forward) primer의 3' end를 target_index에 고정
		# Primer3에서 SEQUENCE_FORCE_LEFT_END는 프라이머의 마지막 염기(3' 말단) 위치를 지정합니다.
		self.update_primer3_seq_args({"SEQUENCE_FORCE_LEFT_END": int(self.target_index)})

		# 3. Primer3 전역 설정 업데이트
		self.update_primer3_global_args(
			{
				"PRIMER_PICK_LEFT_PRIMER": 1,
				"PRIMER_PICK_RIGHT_PRIMER": 1,
				"PRIMER_NUM_RETURN": int(self.n_primers),
				"PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
			}
		)

		# 4. 혼선 방지: 반대 방향(Right) 고정 및 다른 포스 설정 제거
		self.primer3_seq_args.pop("SEQUENCE_FORCE_LEFT_START", None)
		self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_START", None)
		self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_END", None)


	def _configure_primer_reverse_fix(self) -> None:
		super()._configure_primer_common()

		# AS-PCR은 타겟을 피하면 안 되므로 target 회피 해제
		self.primer3_seq_args.pop("SEQUENCE_TARGET", None)

		# ✅ RIGHT primer의 3' end를 target_index에 고정
		# primer3에서 RIGHT의 "start"는 3' end 좌표로 취급되는 것이 일반적
		self.update_primer3_seq_args({"SEQUENCE_FORCE_RIGHT_START": int(self.target_index)})

		self.update_primer3_global_args(
			{
				"PRIMER_PICK_LEFT_PRIMER": 1,
				"PRIMER_PICK_RIGHT_PRIMER": 1,
				"PRIMER_NUM_RETURN": int(self.n_primers),
				"PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
			}
		)

		# 혼선 방지: 반대 force 제거
		self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_END", None)
		self.primer3_seq_args.pop("SEQUENCE_FORCE_LEFT_START", None)
		self.primer3_seq_args.pop("SEQUENCE_FORCE_LEFT_END", None)

	def design(self) -> List[Amplicon]:
		out = []
		out.extend(self._run_mode_and_build("forward_fix"))
		out.extend(self._run_mode_and_build("reverse_fix"))
		self.amplicon_list = out
		return out

	def _run_mode_and_build(self, mode: str) -> List[Amplicon]:
		self.reset()
		print(self.primer3_global_args)
		print(self.primer3_seq_args)
		if mode == "forward_fix": self._configure_primer_forward_fix()
		else: self._configure_primer_reverse_fix()

		res = primer3.bindings.designPrimers(self.primer3_seq_args, self.primer3_global_args)
		print('Results',res)
		return self._build_amplicons_from_result(res or {}, mode=mode)

	def _build_amplicons_from_result(self, res: Dict[str, Any], mode: str) -> List[Amplicon]:
		n_pairs = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))
		if n_pairs == 0: return []
		print(n_pairs)
		amplicons = []
		ref_tpl_raw = self.reference_template_sequence
		alt_tpl_raw = self.template_sequence

		for i in range(n_pairs):
			l_seq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
			r_seq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
			l_s, l_e = _pos_to_binding_strict(res[f"PRIMER_LEFT_{i}"])
			r_s, r_e = _pos_to_binding_strict(res[f"PRIMER_RIGHT_{i}"])

			# 1. Generate All Primer Variants
			# Standard 3' match
			fw_ref = _replace_3prime_base(l_seq, self.ref_allele, "forward")
			fw_alt = _replace_3prime_base(l_seq, self.alt_allele, "forward")
			rv_ref = _replace_3prime_base(r_seq, self.ref_allele, "reverse")
			rv_alt = _replace_3prime_base(r_seq, self.alt_allele, "reverse")
			
			# Specificity Mismatch (pos=3)
			fw_ref_mm = _apply_as_pcr_logic(l_seq, self.ref_allele, "forward", mismatch_pos=3)
			fw_alt_mm = _apply_as_pcr_logic(l_seq, self.alt_allele, "forward", mismatch_pos=3)
			rv_ref_mm = _apply_as_pcr_logic(r_seq, self.ref_allele, "reverse", mismatch_pos=3)
			rv_alt_mm = _apply_as_pcr_logic(r_seq, self.alt_allele, "reverse", mismatch_pos=3)

			# 2. Logic branching by Mode
			if mode == "forward_fix":
				# Forward is Allele-Specific, Right is Common (from Primer3)
				configs = [
					("wt", fw_ref, r_seq, ref_tpl_raw, "ref"),
					("alt", fw_alt, r_seq, alt_tpl_raw, "alt"),
					("wt_mismatch", fw_ref_mm, r_seq, ref_tpl_raw, "ref"),
					("alt_mismatch", fw_alt_mm, r_seq, alt_tpl_raw, "alt"),
				]
			else: # reverse_fix
				# Reverse is Allele-Specific, Left is Common (from Primer3)
				configs = [
					("wt", l_seq, rv_ref, ref_tpl_raw, "ref"),
					("alt", l_seq, rv_alt, alt_tpl_raw, "alt"),
					("wt_mismatch", l_seq, rv_ref_mm, ref_tpl_raw, "ref"),
					("alt_mismatch", l_seq, rv_alt_mm, alt_tpl_raw, "alt"),
				]

			# 3. Create Amplicons
			for label, f_seq, r_seq_final, base_tpl, allele_type in configs:
				# Patch template with chosen primers
				patched_tpl = _patch_template_for_pair(
					base_tpl, left_seq_5to3=f_seq, left_start=l_s, left_end=l_e,
					right_seq_5to3=r_seq_final, right_start=r_s, right_end=r_e
				)

				f_primer = _mk_primer(
					template_sequence=patched_tpl, reference_template_sequence=ref_tpl_raw,
					sequence=f_seq, strand="forward", primer_type="forward",
					target_start_index=self.target_start_index, target_end_index=self.target_end_index,
					binding_start_index=l_s, binding_end_index=l_e
				)
				r_primer = _mk_primer(
					template_sequence=patched_tpl, reference_template_sequence=ref_tpl_raw,
					sequence=r_seq_final, strand="reverse", primer_type="reverse",
					target_start_index=self.target_start_index, target_end_index=self.target_end_index,
					binding_start_index=r_s, binding_end_index=r_e
				)

				amplicons.append(Amplicon(
					template_sequence=patched_tpl,
					reference_template_sequence=ref_tpl_raw,
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					forward_primer=f_primer,
					reverse_primer=r_primer,
					assay=f"as_pcr::{mode}::set{i}::{label}",
					allele=allele_type
				))

		return amplicons