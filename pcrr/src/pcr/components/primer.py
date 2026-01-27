# ... (상단 Import 생략, 기존 Primer 클래스 유지) ...

# [MODIFIED] targets 모듈 임포트
from __future__ import annotations
from .targets import TargetType, Variant, CpG 
from typing import Any, Dict, List, Optional, Literal
from pydantic import BaseModel, Field

import primer3
from Bio.Seq import reverse_complement
from Bio.SeqUtils import gc_fraction

# 패키지 내부 utils 참조
from ..utils import get_start_end_index


Allele = Literal["ref", "alt"]

class Primer(BaseModel):
	def __init__(
		self,
		template_sequence: str,
		sequence: str,
		strand: str,
		primer_type: str,
		target_start_index: int,
		target_end_index: int,
		reference_template_sequence: Optional[str] = None,
		chrom: Optional[str] = None,
		start: Optional[int] = None,
		end: Optional[int] = None,
		binding_start_index: Optional[int] = None,
		binding_end_index: Optional[int] = None,
		tm: Optional[float] = Field(None, description="Melting Temperature of the primer"),
		salt_monovalent_conc: float = 50.0,
		salt_divalent_conc: float = 1.5,
		dntp_conc: float = 0.6,
		dna_conc: float = 50.0,
	) -> None:
		self.template_sequence = template_sequence
		self.reference_template_sequence = reference_template_sequence
		self.sequence = sequence
		self.strand = strand
		self.primer_type = primer_type
		self.target_start_index = target_start_index
		self.target_end_index = target_end_index

		self.length = len(sequence)
		
		# 좌표 설정 로직
		if binding_start_index is not None and binding_end_index is not None:
			self.binding_start_index = binding_start_index
			self.binding_end_index = binding_end_index
		else:
			self.binding_start_index, self.binding_end_index = get_start_end_index(self.template_sequence, self.sequence)
		
		self.chrom = chrom
		self.start = start
		self.end = end

		self.salt_monovalent_conc = salt_monovalent_conc
		self.salt_divalent_conc = salt_divalent_conc
		self.dntp_conc = dntp_conc
		self.dna_conc = dna_conc

		self._calc_basic_properties()
		self._calc_hairpin()
		self._calc_homodimer()

	def _calc_basic_properties(self) -> None:
		self.tm = primer3.calc_tm(
			self.sequence,
			mv_conc=self.salt_monovalent_conc,
			dv_conc=self.salt_divalent_conc,
			dntp_conc=self.dntp_conc,
			dna_conc=self.dna_conc,
		)
		self.gc_percent = gc_fraction(self.sequence, ambiguous="ignore") * 100.0

	def _calc_hairpin(self) -> None:
		result = primer3.calc_hairpin(
			self.sequence,
			mv_conc=self.salt_monovalent_conc,
			dv_conc=self.salt_divalent_conc,
			dntp_conc=self.dntp_conc,
			dna_conc=self.dna_conc,
		)
		self.hairpin = result.structure_found
		self.hairpin_tm = result.tm
		self.hairpin_dg = result.dg / 1000.0
		self.hairpin_dh = result.dh / 1000.0
		self.hairpin_ds = result.ds / 1000.0

	def _calc_homodimer(self) -> None:
		result = primer3.calc_homodimer(
			self.sequence,
			mv_conc=self.salt_monovalent_conc,
			dv_conc=self.salt_divalent_conc,
			dntp_conc=self.dntp_conc,
			dna_conc=self.dna_conc,
		)
		self.homodimer = result.structure_found
		self.homodimer_tm = result.tm
		self.homodimer_dg = result.dg / 1000.0
		self.homodimer_dh = result.dh / 1000.0
		self.homodimer_ds = result.ds / 1000.0

	# -------------------------
	# checks
	# -------------------------
	def check_three_prime_is(self, sequence: str, *, use_reference: bool = True) -> bool:
		ref = self.reference_template_sequence if use_reference else self.template_sequence
		seq_len = len(sequence)

		if self.strand == "forward":
			# end_index가 None일 경우 안전장치 필요 (여기선 있다고 가정)
			end = self.binding_end_index 
			start = max(end - seq_len, 0)
			three_prime_seq = ref[start:end]
		elif self.strand == "reverse":
			start = self.binding_start_index
			end = min(self.binding_start_index + seq_len, len(ref))
			# Reverse strand 3' end corresponds to 5' end of complementary seq in Ref
			three_prime_seq = reverse_complement(ref[start:end])
		else:
			raise ValueError(f"Unknown strand type: {self.strand}")

		return three_prime_seq == sequence
		
	# ... (gc_clamp, count_cpg 등 나머지 메서드 생략 없이 사용하시면 됩니다) ...
		
	def to_dict(self, ignore_attributes: Optional[List[str]] = None) -> Dict[str, Any]:
		if ignore_attributes is None:
			ignore_attributes = [
				"template_sequence",
				"reference_template_sequence",
				"primer_type",
				"chrom",
				"start",
				"end",
				"target_start_index",
				"target_end_index",
			]
		d: Dict[str, Any] = {}
		for k, v in self.__dict__.items():
			if k in ignore_attributes:
				continue
			d[f"{self.primer_type}_{k}"] = v
		return d

class Probe(Primer):
	"""
	SNP/Indel/CpG를 타겟팅하는 Probe.
	Target의 구간(Start~End)을 완전히 포함해야 함.
	"""
	def __init__(self, *args, target: Optional[TargetType] = None, **kwargs) -> None:
		super().__init__(*args, primer_type="probe", **kwargs)
		self.target = target

	def covers_target(self) -> bool:
		"""Probe가 Target의 전체 구간을 포함하는지 확인"""
		if not self.target:
			return False
		
		if self.binding_start_index is None or self.binding_end_index is None:
			return False

		# 조건: Probe 시작점 <= 타겟 시작점  AND  타겟 끝점 <= Probe 끝점
		# (Indel 처리를 위해 전체 포함 여부 확인)
		return (self.binding_start_index <= self.target.index_start) and \
			   (self.target.index_end <= self.binding_end_index)

	def validate_specificity(self) -> bool:
		"""
		Probe Sequence가 타겟(Allele/Methylation) 서열과 일치하는지 확인.
		"""
		if not self.covers_target():
			return False

		# 1. 기대 서열 (Expected String)
		#	Variant(Indel)인 경우 "ATGC" 처럼 여러 글자일 수 있음
		expected_seq = self.target.get_expected_sequence().upper()

		# 2. 실제 Probe 내에서의 타겟 서열 추출
		#	Probe 시작점 기준으로 상대 좌표 계산
		rel_start = self.target.index_start - self.binding_start_index
		rel_end = self.target.index_end - self.binding_start_index
		
		# 범위 체크 (covers_target에서 했지만 안전장치)
		if rel_start < 0 or rel_end > len(self.sequence):
			return False
			
		# Probe Sequence에서 해당 부분 잘라내기
		actual_seq = self.sequence[rel_start:rel_end].upper()

		# 3. 비교
		return actual_seq == expected_seq

	def to_dict(self, ignore_attributes: Optional[List[str]] = None) -> Dict[str, Any]:
		d = super().to_dict(ignore_attributes)
		
		if self.target:
			d["target_type"] = self.target.__class__.__name__
			d["target_chrom"] = self.target.chrom
			
			# Variant인 경우 start/end 기록
			if hasattr(self.target, 'start'):
				d["target_start"] = self.target.start
				d["target_end"] = self.target.end
			else:
				d["target_pos"] = getattr(self.target, 'pos', None)
			
			if isinstance(self.target, Variant):
				d["target_detail"] = f"{self.target.target_allele} ({self.target.ref}>{self.target.alt})"
			elif isinstance(self.target, CpG):
				d["target_detail"] = f"{self.target.methylation_status}"
				
		return d