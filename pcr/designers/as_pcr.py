from __future__ import annotations
import primer3
from typing import Any, Dict, List, Optional, Tuple
from Bio.Seq import Seq

from pcr.designers.base import BasePrimerDesigner
from pcr.components import Amplicon
from pcr.components.primer import Primer

# ----------------------------------------------------------------------
# Helper Functions for AS-PCR Logic
# ----------------------------------------------------------------------

def _replace_3prime_base(primer_5to3: str, base: str, strand: str = "forward") -> str:
	"""Replaces the 3' end of the primer with the target allele."""
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
	"""Applies the 3' allele and an internal artificial mismatch."""
	mismatch_map = {
		"strong": {'A': 'G', 'T': 'C', 'G': 'A', 'C': 'T'},
		"medium": {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'},
		"weak":   {'A': 'C', 'T': 'G', 'G': 'T', 'C': 'A'},
	}
	primer_seq = _replace_3prime_base(primer_5to3, target_base, strand)
	s = list(primer_seq)
		
	idx = -int(mismatch_pos) # Position from 3' end
	if abs(idx) <= len(s):
		original = s[idx]
		new_base = mismatch_map.get(intensity.lower(), mismatch_map["strong"]).get(original, "A")
		if new_base == original:
			for cand in ("A", "C", "G", "T"):
				if cand != original:
					new_base = cand
					break
		s[idx] = new_base
	return "".join(s)

def _patch_template_for_pair(
	template_5to3: str,
	reference_template_5to3: str,
	*,
	left_seq_5to3: str,
	left_start: int,
	left_end: int,
	right_seq_5to3: str,
	right_start: int,
	right_end: int,
) -> Tuple[str, List[Dict[str, Any]]]:
	"""Overwrites template with primers and records diffs from reference."""
	t_list = list(template_5to3.upper())
	ref_list = list(reference_template_5to3.upper())
		
	left_seq = left_seq_5to3.upper()
	right_rc = str(Seq(right_seq_5to3.upper()).reverse_complement())

	# Apply primer sequences to the template array
	t_list[left_start : left_end + 1] = list(left_seq)
	t_list[right_start : right_end + 1] = list(right_rc)
	patched_template = "".join(t_list)

	# Compare base-by-base to identify SNP vs Artificial Mismatches
	diffs = []
	for idx in range(min(len(t_list), len(ref_list))):
		if t_list[idx] != ref_list[idx]:
			source = "GENOMIC_VARIANT" # Default: the SNP itself
			if left_start <= idx <= left_end:
				source = "FORWARD_PRIMER_PATCH"
			elif right_start <= idx <= right_end:
				source = "REVERSE_PRIMER_PATCH"
			
			diffs.append({
				"index": idx,
				"ref": ref_list[idx],
				"patched": t_list[idx],
				"source": source
			})
	return patched_template, diffs

# ----------------------------------------------------------------------
# Main AS-PCR Designer Class
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

	def _configure_primer_forward_fix(self) -> None:
		"""Forces the Forward primer to end (3') at the target SNP."""
		super()._configure_primer_common()
		self.primer3_seq_args.pop("SEQUENCE_TARGET", None)
		# 3' end of the Forward primer = target_index
		self.update_primer3_seq_args({"SEQUENCE_FORCE_LEFT_END": self.target_index})
		self.update_primer3_global_args({
			"PRIMER_PICK_LEFT_PRIMER": 1,
			"PRIMER_PICK_RIGHT_PRIMER": 1,
			"PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
		})

	def _configure_primer_reverse_fix(self) -> None:
		"""Forces the Reverse primer to start (3') at the target SNP."""
		super()._configure_primer_common()
		self.primer3_seq_args.pop("SEQUENCE_TARGET", None)
		# In Primer3, 'FORCE_RIGHT_END' defines the 3' end of the Reverse primer
		self.update_primer3_seq_args({"SEQUENCE_FORCE_RIGHT_END": self.target_index})
		self.update_primer3_global_args({
			"PRIMER_PICK_LEFT_PRIMER": 1,
			"PRIMER_PICK_RIGHT_PRIMER": 1,
			"PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
		})

	def design(self) -> List[Amplicon]:
		results = []
		# Run design for both orientations
		results.extend(self._run_mode_and_build("forward_fix"))
		results.extend(self._run_mode_and_build("reverse_fix"))
		self.amplicon_list = results
		return results

	def _run_mode_and_build(self, mode: str) -> List[Amplicon]:
		self.reset()
		if mode == "forward_fix":
			self._configure_primer_forward_fix()
		else:
			self._configure_primer_reverse_fix()

		res = primer3.bindings.designPrimers(self.primer3_seq_args, self.primer3_global_args)
		return self._build_amplicons_from_result(res or {}, mode=mode)

	def _build_amplicons_from_result(self, res: Dict[str, Any], mode: str) -> List[Amplicon]:
		num_returned = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
		amplicons = []
		
		for i in range(num_returned):
			# 1. Extract raw Primer3 data
			l_seq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
			r_seq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
			l_pos = res[f"PRIMER_LEFT_{i}"] # (start_index, length)
			r_pos = res[f"PRIMER_RIGHT_{i}"] # (3_prime_index, length)
			
			# 2. COORDINATE FIX
			# Forward: Starts at index, ends at index + len - 1
			l_s, l_e = l_pos[0], l_pos[0] + l_pos[1] - 1
			
			# Reverse: r_pos[0] IS the 3' end (the lowest index on the sense strand).
			# The binding site spans from r_pos[0] to r_pos[0] + length - 1.
			r_len = r_pos[1]
			r_binding_start = r_pos[0]
			r_binding_end = r_pos[0] + r_len - 1

			# 3. Define configurations based on the fixed primer orientation
			if mode == "forward_fix":
				# In this mode, only the Forward primer 3' end is at target_index
				configs = [
					("wt", _replace_3prime_base(l_seq, self.ref_allele, "forward"), r_seq, self.reference_template_sequence, "ref"),
					("alt", _replace_3prime_base(l_seq, self.alt_allele, "forward"), r_seq, self.template_sequence, "alt"),
					("wt_mm", _apply_as_pcr_logic(l_seq, self.ref_allele, "forward"), r_seq, self.reference_template_sequence, "ref"),
					("alt_mm", _apply_as_pcr_logic(l_seq, self.alt_allele, "forward"), r_seq, self.template_sequence, "alt")
				]
			else: # reverse_fix
				# In this mode, only the Reverse primer 3' end is at target_index
				configs = [
					("wt", l_seq, _replace_3prime_base(r_seq, self.ref_allele, "reverse"), self.reference_template_sequence, "ref"),
					("alt", l_seq, _replace_3prime_base(r_seq, self.alt_allele, "reverse"), self.template_sequence, "alt"),
					("wt_mm", l_seq, _apply_as_pcr_logic(r_seq, self.ref_allele, "reverse"), self.reference_template_sequence, "ref"),
					("alt_mm", l_seq, _apply_as_pcr_logic(r_seq, self.alt_allele, "reverse"), self.template_sequence, "alt")
				]

			for label, f_seq, r_seq_final, base_tpl, allele_type in configs:
				# 4. Patch template and identify differences
				# The Reverse primer sequence (r_seq_final) is 5'->3'.
				# _patch_template_for_pair will RC it and place it at [r_binding_start : r_binding_end + 1]
				patched_tpl, diffs = _patch_template_for_pair(
					base_tpl, self.reference_template_sequence,
					left_seq_5to3=f_seq, left_start=l_s, left_end=l_e,
					right_seq_5to3=r_seq_final, right_start=r_binding_start, right_end=r_binding_end
				)

				# 5. Correctly initialize Primer objects
				f_primer = Primer(
					template_sequence=patched_tpl,
					primer_type="forward",
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					sequence=f_seq,
					strand="forward",
					binding_start_index=l_s,
					binding_end_index=l_e
				)

				r_primer = Primer(
					template_sequence=patched_tpl,
					primer_type="reverse",
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					sequence=r_seq_final,
					strand="reverse",
					binding_start_index=r_binding_start,
					binding_end_index=r_binding_end
				)

				# 6. Build and collect the Amplicon
				amp = Amplicon(
					template_sequence=patched_tpl,
					reference_template_sequence=self.reference_template_sequence,
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					forward_primer=f_primer,
					reverse_primer=r_primer,
					assay=f"AS-PCR::{mode}::{label}",
					allele=allele_type
				)
				amplicons.append(amp)
				
		return amplicons