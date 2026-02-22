from __future__ import annotations

from typing import Any, Dict, List, Optional, Literal
from dataclasses import dataclass, asdict

@dataclass
class GenomicRegion:
	"""
	유전체 상의 위치 정보를 담는 객체.
	Primer나 Probe, Amplicon이 이 객체를 멤버로 가짐.
	"""
	chrom: str
	start: int  # 0-based start
	end: int	# 0-based end (exclusive)
	strand: str # "+" or "-"
	sequence: str = ''
	base: int = 0
	
	def __len__(self) -> int:
		return self.end - self.start

	def to_dict(self, prefix: str = "") -> Dict[str, Any]:
		"""Flatten하여 반환 (예: forward_primer_chrom)"""
		p = f"{prefix}_" if prefix else ""
		region_id = f"{self.chrom}:{self.start}-{self.end}"
		return {
			f"{p}chrom": self.chrom,
			f"{p}start": self.start,
			f"{p}end": self.end,
			f"{p}strand": self.strand,
			f"{p}base": self.base,
			f"{p}id": region_id
		}

	def get_id(self):
		return f"{self.chrom}:{self.start}-{self.end}"

@dataclass
class SequenceChange:
	"""
	Template(Target)과 Reference 사이의 염기 서열 차이를 기록.
	"""
	position: int	  # Amplicon 내 상대 위치 (0-based)
	ref_base: str	  # Reference 염기
	alt_base: str	  # Template 염기
	region_type: str   # 'forward_primer', 'internal', etc.
	on_target: bool	# 타겟 여부
		
	# ✅ 이름 명확화: 단순 변환이 아니라 'Bisulfite Conversion' 패턴인지
	is_bisulfite_conversion: bool 

	def to_dict(self) -> Dict[str, Any]:
		return {
			"pos": self.position,
			"region": self.region_type,
			"ref": self.ref_base,
			"alt": self.alt_base,
			"change": f"{self.ref_base}>{self.alt_base}",
			"on_target": self.on_target,
			"is_conversion": self.is_bisulfite_conversion
		}