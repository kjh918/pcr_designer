from __future__ import annotations

from typing import Any, Dict, Optional, Literal

from pcr.components.primer import Primer
from pcr.utils import get_start_end_index
import primer3
from Bio.SeqUtils import gc_fraction
#from pcr.seq.fetch import 


class Amplicon:
	def __init__(
		self,
		template_sequence: str,
		target_start_index: int,
		target_end_index: int,
		reference_template_sequence: Optional[str] = None,
		chrom: Optional[str] = None,
		start: Optional[int] = None,
		end: Optional[int] = None,
		forward_primer: Optional[Primer] = None,
		reverse_primer: Optional[Primer] = None,
		probe: Optional[Primer] = None,
		assay: str = "generic",
		allele: Allele = "ref",
	) -> None:
		self.template_sequence = template_sequence
		self.reference_template_sequence = reference_template_sequence # or template_sequence

		self.target_start_index = target_start_index
		self.target_end_index = target_end_index

		self.chrom = chrom
		self.start = start
		self.end = end

		self.forward_primer = forward_primer
		self.reverse_primer = reverse_primer
		self.probe = probe

		self.assay = assay
		self.allele = allele

		if self.forward_primer is not None:
			if self.forward_primer.binding_start_index is not None and self.forward_primer.binding_end_index is not None:
				self.forward_start_index, self.forward_end_index = self.forward_primer.binding_start_index, self.forward_primer.binding_end_index
			else:
				self.forward_start_index, self.forward_end_index = get_start_end_index(
					self.forward_primer.template_sequence, self.forward_primer.sequence
				)
		if self.reverse_primer is not None:
			if self.reverse_primer.binding_start_index is not None and self.reverse_primer.binding_end_index is not None:
				self.reverse_start_index, self.reverse_end_index = self.reverse_primer.binding_start_index, self.reverse_primer.binding_end_index
			else:
				self.reverse_start_index, self.reverse_end_index = get_start_end_index(
					self.reverse_primer.template_sequence, self.reverse_primer.sequence
				)
		if self.probe is not None:		
			if self.probe.binding_start_index is not None and self.probe.binding_end_index is not None:
				self.probe_start_index, self.probe_end_index = self.probe.binding_start_index, self.probe.binding_end_index
			else:
				self.probe_start_index, self.probe_end_index = get_start_end_index(
					self.probe.template_sequence, self.probe.sequence
				)
			
		# Amplicon sequence (template 기준)
		self.amplicon_sequence: Optional[str] = self._calc_amplicon_sequence()

		# Amplicon metrics (template 기준)
		self.amplicon_gc: Optional[float] = None
		self.amplicon_tm: Optional[float] = None

		# Reference 기준 (template/reference가 다를 때만 의미있음)
		self.reference_amplicon_sequence: Optional[str] = None
		self.reference_amplicon_gc: Optional[float] = None
		self.reference_amplicon_tm: Optional[float] = None

		self._calc_amplicon_metrics()

	# -------------------------
	# helpers
	# -------------------------
	@staticmethod
	def _calc_tm_primer3(seq: str) -> float:
		s = (seq or "").upper()
		if not s:
			raise ValueError("Empty sequence for calcTm")
		return float(primer3.calc_tm(s))

	@staticmethod
	def _calc_gc_primer3(seq: str) -> float:
		s = (seq or "").upper()
		return gc_fraction(s, ambiguous="ignore") * 100.0
	# -------------------------
	# core
	# -------------------------
	def _calc_amplicon_sequence(self) -> Optional[str]:
		if self.forward_primer is None or self.reverse_primer is None:
			return None	
		return self.template_sequence[self.forward_start_index : self.reverse_end_index + 1]

	def _calc_reference_amplicon_sequence(self) -> Optional[str]:
		if self.forward_primer is None or self.reverse_primer is None:
			return None
		ref = self.reference_template_sequence
		if not ref:
			return None

		# template에서 계산된 primer 좌표를 reference에도 그대로 적용
		if self.reverse_end_index >= len(ref) or self.forward_start_index < 0:
			# indel 등으로 길이가 달라져 인덱스가 깨진 경우 방어
			return None
		return ref[self.forward_start_index : self.reverse_end_index + 1]

	def _calc_amplicon_metrics(self) -> None:
		# template 기준
		if self.amplicon_sequence:
			self.amplicon_gc = self._calc_gc_primer3(self.amplicon_sequence)
			self.amplicon_tm = self._calc_tm_primer3(self.amplicon_sequence)

		# reference 고려 (reference != template 인 경우)
		if self.reference_template_sequence != self.template_sequence:
			self.reference_amplicon_sequence = self._calc_reference_amplicon_sequence()
			if self.reference_amplicon_sequence:
				self.reference_amplicon_gc = self._calc_gc_primer3(self.reference_amplicon_sequence)
				self.reference_amplicon_tm = self._calc_tm_primer3(self.reference_amplicon_sequence)

	def to_dict(self) -> Dict[str, Any]:
		d: Dict[str, Any] = {
			"reference_template_sequence": self.reference_template_sequence,
			"template_sequence": self.template_sequence,
			"target_start_index": self.target_start_index,
			"target_end_index": self.target_end_index,
			"assay": self.assay,
			"allele": self.allele,
		}

		if self.amplicon_sequence is not None:
			d["amplicon_sequence"] = self.amplicon_sequence
			d["amplicon_length"] = len(self.amplicon_sequence)
		else:
			d["amplicon_sequence"] = None
			d["amplicon_length"] = None

		# ✅ Amplicon metrics (template)
		d["amplicon_gc"] = self.amplicon_gc
		d["amplicon_tm"] = self.amplicon_tm

		# ✅ Reference-aware metrics
		# (reference가 같으면 None으로 두거나, 같을 때도 채우고 싶으면 조건을 빼면 됨)
		d["reference_amplicon_sequence"] = self.reference_amplicon_sequence
		d["reference_amplicon_gc"] = self.reference_amplicon_gc
		d["reference_amplicon_tm"] = self.reference_amplicon_tm

		if self.forward_primer is not None:
			d.update(self.forward_primer.to_dict())
		if self.reverse_primer is not None:
			d.update(self.reverse_primer.to_dict())
		if self.probe is not None:
			d.update(self.probe.to_dict())
		return d
