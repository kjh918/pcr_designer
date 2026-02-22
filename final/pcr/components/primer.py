from __future__ import annotations

from typing import Any, Dict, List, Optional, Literal
from dataclasses import dataclass, field
import primer3

## ##
from .region import GenomicRegion


ROLE_MAP = {"LEFT": "forward", "RIGHT": "reverse", "INTERNAL": "probe"}
# ==============================================================================
# 1. Primer
# ==============================================================================
class Primer:
	def __init__(
		self,
		sequence: str, 
		role: str,
		length: int,
		penalty: float,
		
		tm: float,
		gc_percent: float,
		hairpin_tm: float,
		hairpin_dg: float,
		homodimer_tm: float,
		homodimer_dg: float,

		cpg_count: int,			   
		start_index: int,	 # 0-based template start
		end_index: int,	   # 0-based template end (exclusive)
		region: Optional[GenomicRegion] = None,
		calc_args: Dict[str, float] = {'mv_conc': 50, 'dv_conc': 1.5, 'dntp_conc': 0.6, 'dna_conc': 50}
	) -> None:
		self.sequence = sequence
		self.role = role
		self.length = length
		self.tm = tm
		self.hairpin_tm = hairpin_tm
		self.hairpin_dg = hairpin_dg
		self.homodimer_tm = homodimer_tm
		self.homodimer_dg = homodimer_dg
		self.gc_percent = gc_percent
		self.penalty = penalty
		self.cpg_count = cpg_count 
		self.start_index = start_index
		self.end_index = end_index
		self.region = region
		self.calc_args = calc_args

	@classmethod
	def from_primer3(
		cls, 
		result: Dict[str, Any], 
		rank: int, 
		role_key: Literal["LEFT", "RIGHT", "INTERNAL"],
		template_region: Optional[GenomicRegion] = None
	) -> Optional['Primer']:
		
		prefix = f"PRIMER_{role_key}_{rank}"
		seq = result.get(f"{prefix}_SEQUENCE")
		if not seq: return None

		# 1. 기본 정보 추출
		role_name = ROLE_MAP.get(role_key, "unknown")
		length = len(seq)
		tm = float(result.get(f"{prefix}_TM", 0.0))
		gc = float(result.get(f"{prefix}_GC_PERCENT", 0.0))
		penalty = float(result.get(f"{prefix}_PENALTY", 0.0))
		cpg_count = seq.count("CG")

		# 2. 2차 구조 정밀 계산 (Hairpin & Homodimer)
		calc_args = {'mv_conc': 50, 'dv_conc': 1.5, 'dntp_conc': 0.6, 'dna_conc': 50}
		# Hairpin Calculation
		hp = primer3.calc_hairpin(seq, **calc_args)
		hp_tm = hp.tm
		hp_dg = hp.dg / 1000.0 if hp.structure_found else 0.0

		# Homodimer Calculation
		hd = primer3.calc_homodimer(seq, **calc_args)
		hd_tm = hd.tm
		hd_dg = hd.dg / 1000.0 if hd.structure_found else 0.0

		# 3. 좌표 계산
		info = result.get(prefix)
		p3_index = info[0] if info else 0

		if role_key == "RIGHT":
			start_index = p3_index - length + 1
		else:
			start_index = p3_index
		end_index = start_index + length

		# 4. Genomic Region 매핑
		region = None
		if template_region:
			if template_region.strand != "-":
				g_start = template_region.start + start_index
				g_end = template_region.start + end_index
			else:
				g_start = template_region.end - end_index
				g_end = template_region.end - start_index
			
			region = GenomicRegion(template_region.chrom, g_start, g_end, template_region.strand)
		return cls(
			sequence=seq, role=role_name, length=length, tm=tm, gc_percent=gc, penalty=penalty, cpg_count=cpg_count,
			hairpin_tm=hp_tm, hairpin_dg=hp_dg,
			homodimer_tm=hd_tm, homodimer_dg=hd_dg,
			start_index=start_index, end_index=end_index, region=region
		)

	def to_dict(self) -> Dict[str, Any]:
		prefix = self.role
		data = {}
		data[f"{prefix}_sequence"] = self.sequence
		data[f"{prefix}_length"] = self.length
		data[f"{prefix}_tm"] = self.tm
		data[f"{prefix}_gc"] = self.gc_percent
		data[f"{prefix}_cpg_count"] = self.cpg_count
		data[f"{prefix}_penalty"] = self.penalty
		
		# ✅ 2차 구조 정보 출력
		data[f"{prefix}_hairpin_tm"] = self.hairpin_tm
		data[f"{prefix}_hairpin_dg"] = self.hairpin_dg
		data[f"{prefix}_homodimer_tm"] = self.homodimer_tm
		data[f"{prefix}_homodimer_dg"] = self.homodimer_dg
		
		if self.region:
			data.update(self.region.to_dict(prefix=prefix))
		return data

@dataclass
class TargetAnnotation:
	id: str
	type: str
	status: str
	region: Optional[GenomicRegion] = None
	metadata: Dict[str, Any] = field(default_factory=dict) 

	def to_dict(self, prefix: str = "probe_target") -> Dict[str, Any]:
		data = {}
		data[f"{prefix}_id"] = self.id
		data[f"{prefix}_type"] = self.type
		data[f"{prefix}_status"] = self.status
		
		if self.region:
			data.update(self.region.to_dict(prefix=prefix))

		for k, v in self.metadata.items():
			data[f"{prefix}_{k}"] = v
		return data

# ==============================================================================
# 2. Probe
# ==============================================================================
class Probe(Primer):
	def __init__(self, target: Optional[TargetAnnotation] = None, **kwargs):
		super().__init__(**kwargs)
		self.target = target

	@classmethod
	def from_primer3(
		cls, 
		result: Dict[str, Any], 
		rank: int, 
		template_region: Optional[GenomicRegion] = None,
		target_id: str = "Unknown",
		target_type: str = "Generic",
		target_status: str = "NA",
		target_region: Optional[GenomicRegion] = None,
		**target_metadata
	) -> Optional['Probe']:
		
		# 1. 부모 메서드를 통해 객체 생성 (cls가 Probe이므로 Probe 인스턴스가 반환됨)
		p = super().from_primer3(result, rank, "INTERNAL", template_region)
		
		if p:
			# 2. TargetAnnotation 생성
			annotation = TargetAnnotation(
				id=target_id, 
				type=target_type, 
				status=target_status,
				region=target_region, 
				metadata=target_metadata
			)
			
			# 3. [수정됨] 객체를 새로 만들지 않고, 생성된 객체의 target 속성만 설정
			p.target = annotation
			return p
			
		return None

	def to_dict(self) -> Dict[str, Any]:
		data = super().to_dict()
		if self.target:
			data.update(self.target.to_dict())
		else:
			data["probe_target_id"] = None
			data["probe_target_type"] = None
		return data