from pydantic import BaseModel, Field, computed_field
from typing import Literal, Dict, Any, Optional

# 개별 Checker의 결과 (Pass 여부 + 데이터 + 에러메시지)
class QCResult(BaseModel):
	passed: bool
	data: Dict[str, Any] = {}
	error: Optional[str] = None

# Amplicon에 붙을 최종 QC 상태 객체
class AmpliconQCStatus(BaseModel):
	is_pass: bool	   # 최종 통과 여부 (AND 조건)
	thermo: QCResult	# 물성(Hairpin/Dimer) 결과
	blast: QCResult	 # BLAST 결과
	ispcr: QCResult	 # isPCR 결과

	def to_summary(self) -> Dict[str, Any]:
		"""Summary용 간단 딕셔너리 반환"""
		return {
			"Final": "PASS" if self.is_pass else "FAIL",
			"Thermo": "OK" if self.thermo.passed else "No",
			"BLAST": "OK" if self.blast.passed else "No",
			"isPCR": "OK" if self.ispcr.passed else "No",
			"Off_Targets(BLAST)": self.blast.data.get("off_target_count", 0),
			"Off_Targets(isPCR)": self.ispcr.data.get("count", 0) - 1 if self.ispcr.passed else "N/A"
		}

class BlastHit(BaseModel):
	"""BLAST 결과 한 줄에 대한 객체"""
	qseqid: str
	sseqid: str  # Chromosome
	pident: float
	length: int
	qstart: int
	qend: int
	sstart: int
	send: int
		
	@computed_field
	def strand(self) -> Literal["+", "-"]:
		return "+" if self.sstart <= self.send else "-"

	@computed_field
	def genomic_start(self) -> int:
		"""Genome 상의 실제 시작 좌표 (항상 작은 값)"""
		return min(self.sstart, self.send)

	@computed_field
	def genomic_end(self) -> int:
		"""Genome 상의 실제 끝 좌표 (항상 큰 값)"""
		return max(self.sstart, self.send)

class OffTargetAmplicon(BaseModel):
	"""Off-target으로 생성될 가능성이 있는 Amplicon 정보"""
	chrom: str
	start: int
	end: int
	product_size: int
	fwd_hit: BlastHit
	rev_hit: BlastHit
		
	@computed_field
	def is_target(self) -> bool:
		"""이 앰플리콘이 우리가 원하는 타겟(Target)인지 판별 (로직은 필요에 따라 확장)"""
		# 보통 pident가 100이고, 쿼리 전체가 매칭되면 타겟으로 간주
		return (self.fwd_hit.pident == 100.0 and self.rev_hit.pident == 100.0)