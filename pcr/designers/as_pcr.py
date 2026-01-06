from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Optional, Tuple

import primer3
from Bio.Seq import Seq

from pcr.designers.base import BasePrimerDesigner
from pcr.components import Amplicon
from pcr.components.primer import Primer


# =============================================================================
# Types / Constants
# =============================================================================

TemplateType = Literal["wt", "alt", "wt_mm", "alt_mm"]
FixedPrime = Literal["forward", "reverse"]
Strand = Literal["forward", "reverse"]


# =============================================================================
# Helpers (coordinate-safe, template-slice based)
# =============================================================================

def _rc(seq: str) -> str:
	return str(Seq(seq).reverse_complement())


def parse_primer_from_template(
	template_5to3: str,
	start: int,
	end: int,
	strand: Strand,
) -> str:
	"""
	Parse primer sequence from template by coordinates.

	template_5to3: sense strand 5'->3'
	start/end: 0-based inclusive
	strand:
	  - "forward": primer sequence is the slice itself
	  - "reverse": primer sequence is reverse-complement of slice
	"""
	if start < 0 or end < 0 or end < start or end >= len(template_5to3):
		raise ValueError(f"Invalid slice: start={start}, end={end}, len={len(template_5to3)}")

	seg = template_5to3[start:end + 1].upper()
	if strand == "forward":
		return seg
	return seg


def primer3_left_pos_to_span(pos: Any) -> Tuple[int, int]:
	"""
	primer3 LEFT position is typically: (start, length)
	-> span: (start, end inclusive)
	"""
	if not pos or len(pos) != 2:
		raise ValueError(f"Invalid LEFT position: {pos}")
	start, length = int(pos[0]), int(pos[1])
	end = start + length - 1
	return start, end


def primer3_right_pos_to_span(pos: Any) -> Tuple[int, int]:
	"""
	primer3 RIGHT position in primer3-py is commonly: (3' end index, length)
	In your convention (as discussed), template binding span is:
	  start = three_prime_index
	  end   = three_prime_index + length - 1
	"""
	if not pos or len(pos) != 2:
		raise ValueError(f"Invalid RIGHT position: {pos}")
	three_prime, length = int(pos[0]), int(pos[1])
	start = three_prime
	end = three_prime + length - 1
	return start, end


def validate_fixed_prime_anchor(
	*,
	fixed_prime: FixedPrime,
	left_span: Tuple[int, int],
	right_span: Tuple[int, int],
	target_index: int,
) -> bool:
	"""
	Enforce:
	  - forward fixed: LEFT 3' end must be target_index -> left_end == target_index
	  - reverse fixed: RIGHT 3' end must be target_index -> right_start == target_index
	"""
	l_s, l_e = left_span
	r_s, r_e = right_span

	if fixed_prime == "forward":
		return l_e == target_index
	else:
		return r_s == target_index


# =============================================================================
# Grouping (wt/alt/wt_mm/alt_mm as a set)
# =============================================================================

@dataclass
class AspcrSet:
	set_id: str				  # e.g. "set0"
	fixed_prime: FixedPrime	  # "forward" or "reverse"
	left_span: Tuple[int, int]   # (start,end)
	right_span: Tuple[int, int]  # (start,end)
	amplicons: Dict[TemplateType, Amplicon]  # keys: wt, alt, wt_mm, alt_mm


# =============================================================================
# AsPcrDesigner (single run of primer3 + parse primers from templates)
# =============================================================================

class AsPcrDesigner(BasePrimerDesigner):
	"""
	Concept:
	  - Run primer3 ONLY on reference_template_sequence ("wt")
	  - Use primer3 coordinates to parse primers from:
		wt / alt / wt_mm / alt_mm templates
	  - Build a 4-amplicon set per primer pair (if anchor constraint satisfied)
	"""

	def __init__(
		self,
		*,
		reference_template_sequence: str,   # wt
		alt_template_sequence: str,		 # alt
		ref_mm_template_sequence: str,	  # wt_mm
		alt_mm_template_sequence: str,	  # alt_mm

		target_start_index: int,
		target_end_index: int,
		target_index: int,

		ref_allele: str,
		alt_allele: str,

		chrom: str,
		start: int,
		end: int,

		mismatch_pos: int = 3,			  # kept for metadata; templates already include it
		fixed_prime: FixedPrime = "forward",# forward / reverse

		# BasePrimerDesigner kwargs
		**kwargs: Any,
	) -> None:
		self.templates: Dict[TemplateType, str] = {
			"wt": reference_template_sequence,
			"alt": alt_template_sequence,
			"wt_mm": ref_mm_template_sequence,
			"alt_mm": alt_mm_template_sequence,
		}

		self.target_index = int(target_index)
		self.ref_allele = ref_allele.upper()
		self.alt_allele = alt_allele.upper()

		self.chrom = chrom
		self.start = int(start)
		self.end = int(end)

		self.mismatch_pos = int(mismatch_pos)
		self.fixed_prime: FixedPrime = fixed_prime

		# IMPORTANT: primer3 must run on WT reference template only
		super().__init__(
			template_sequence=reference_template_sequence,
			reference_template_sequence=reference_template_sequence,
			target_start_index=int(target_start_index),
			target_end_index=int(target_end_index),
			**kwargs,
		)

	# ------------------------------------------------------------------
	# primer3 configuration
	# ------------------------------------------------------------------
	def _configure_force_anchor(self) -> None:
		"""
		Enforce 3' anchor at target_index.
		- forward fixed => force LEFT_END = target_index
		- reverse fixed => force RIGHT_END = target_index (commonly 3' end)
		  But in your later logic, you validated RIGHT_START == target_index (3' index).
		  So here we must align primer3 forcing with your interpretation.

		For primer3:
		  - SEQUENCE_FORCE_LEFT_END forces left primer 3' end (OK)
		  - For right primer, depending on primer3 build, use:
			  SEQUENCE_FORCE_RIGHT_START  (if primer3-py treats RIGHT pos[0] as 3' end)
			or
			  SEQUENCE_FORCE_RIGHT_END	(if primer3 expects end coordinate)
		You have been using pos[0] as 3' end and binding span as [pos[0] : pos[0]+len-1],
		therefore we should force RIGHT_START to be target_index.

		If your primer3 build only supports RIGHT_END, switch accordingly.
		"""
		super()._configure_primer_common()

		# AS-PCR must not avoid target
		self.primer3_seq_args.pop("SEQUENCE_TARGET", None)

		if self.fixed_prime == "forward":
			self.update_primer3_seq_args({"SEQUENCE_FORCE_LEFT_END": self.target_index})
			# clear possible conflicting right constraints
			self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_START", None)
			self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_END", None)
		else:
			# IMPORTANT: aligns with our RIGHT span convention (pos[0] == 3' end)
			self.update_primer3_seq_args({"SEQUENCE_FORCE_RIGHT_START": self.target_index})
			self.primer3_seq_args.pop("SEQUENCE_FORCE_LEFT_END", None)
			self.primer3_seq_args.pop("SEQUENCE_FORCE_LEFT_START", None)

		# Ensure both primers are picked
		self.update_primer3_global_args(
			{
				"PRIMER_PICK_LEFT_PRIMER": 1,
				"PRIMER_PICK_RIGHT_PRIMER": 1,
				"PRIMER_NUM_RETURN": int(self.n_primers),
				"PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
				"PRIMER_EXPLAIN_FLAG": 1,
			}
		)

	# ------------------------------------------------------------------
	# main design entry
	# ------------------------------------------------------------------
	def design_sets(self) -> List[AspcrSet]:
		"""
		Returns list of AspcrSet (each contains wt/alt/wt_mm/alt_mm amplicons).
		Also populates self.amplicon_list as a flat list (all amplicons).
		"""
		self.reset()
		self._configure_force_anchor()

		res = primer3.bindings.designPrimers(self.primer3_seq_args, self.primer3_global_args) or {}
		n_pairs = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))

		if n_pairs == 0:
			# keep explain for debugging
			print("[Primer3 Explain]",
				  res.get("PRIMER_LEFT_EXPLAIN"),
				  res.get("PRIMER_RIGHT_EXPLAIN"),
				  res.get("PRIMER_PAIR_EXPLAIN"))
			self.amplicon_list = []
			return []

		sets: List[AspcrSet] = []
		flat: List[Amplicon] = []
		for i in range(n_pairs):
			left_pos = res.get(f"PRIMER_LEFT_{i}")
			right_pos = res.get(f"PRIMER_RIGHT_{i}")
			if left_pos is None or right_pos is None:
				continue

			left_span = primer3_left_pos_to_span(left_pos)
			right_span = primer3_right_pos_to_span(right_pos)

			# (optional) anchor check
			if not validate_fixed_prime_anchor(
				fixed_prime=self.fixed_prime,
				left_span=left_span,
				right_span=right_span,
				target_index=self.target_index,
			):
				continue

			set_id = f"set{i}"
			amplicons_by_type = {}

			for ttype in ("wt", "alt", "wt_mm", "alt_mm"):
				tpl = self.templates[ttype]

				# ✅ primer sequence는 primer3 output 말고 "index로 template에서" 가져온다
				f_seq = parse_primer_from_template(tpl, left_span[0], left_span[1], "forward")
				r_seq = parse_primer_from_template(tpl, right_span[0], right_span[1], "reverse")

				f_primer = Primer(
					template_sequence=tpl,
					reference_template_sequence=self.reference_template_sequence,  # genome ref 쓰고 싶으면 그걸로
					sequence=f_seq,
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					strand="forward",
					primer_type="forward",
					binding_start_index=left_span[0],
					binding_end_index=left_span[1],
				)
				r_primer = Primer(
					template_sequence=tpl,
					reference_template_sequence=self.reference_template_sequence,
					sequence=r_seq,
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					strand="reverse",
					primer_type="reverse",
					binding_start_index=right_span[0],
					binding_end_index=right_span[1],
				)

				amp = Amplicon(
					template_sequence=tpl,
					reference_template_sequence=self.reference_template_sequence,
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					forward_primer=f_primer,
					reverse_primer=r_primer,
					assay=f"as_pcr::{self.fixed_prime}::{set_id}::{ttype}",
					allele=ttype,
				)
				amplicons_by_type[ttype] = amp
				flat.append(amp)


		self.amplicon_list = flat
		return sets

	def design(self) -> List[Amplicon]:
		"""
		Compatibility: returns flat list (like other designers).
		"""
		self.design_sets()
		return self.amplicon_list
