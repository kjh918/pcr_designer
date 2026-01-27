from pydantic import BaseModel, Field
from typing import Dict, Any, Optional

# 공통 설정을 위한 기본 클래스
class BaseDesignKwargs(BaseModel):
	opt_length: int = 20
	min_length: int = 18
	max_length: int = 27
	opt_tm: float = 60.0
	min_tm: float = 55.0
	max_tm: float = 65.0
	opt_gc: float = 45.0
	min_gc: float = 25.0
	max_gc: float = 85.0
	n_primers: int = 100

	def to_primer3_args(self, prefix: str = "PRIMER_") -> Dict[str, Any]:
		"""YAML의 friendly name을 Primer3 옵션명으로 변환"""
		return {
			f"{prefix}OPT_SIZE": self.opt_length,
			f"{prefix}MIN_SIZE": self.min_length,
			f"{prefix}MAX_SIZE": self.max_length,
			f"{prefix}OPT_TM": self.opt_tm,
			f"{prefix}MIN_TM": self.min_tm,
			f"{prefix}MAX_TM": self.max_tm,
			f"{prefix}OPT_GC_PERCENT": self.opt_gc,
			f"{prefix}MIN_GC": self.min_gc,
			f"{prefix}MAX_GC": self.max_gc,
			f"{prefix}NUM_RETURN": self.n_primers,
		}

class PrimerKwargs(BaseDesignKwargs):
	min_amplicon_length: int = 80
	max_amplicon_length: int = 120
	n_primers: int = 100
	primer3_global_args: Dict[str, Any] = Field(default_factory=dict)

	def to_global_args(self) -> Dict[str, Any]:
		"""Primer 설계를 위한 전체 Global Args 생성"""
		args = self.to_primer3_args(prefix="PRIMER_")
		args.update({
			"PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
			"PRIMER_NUM_RETURN": self.n_primers,
		})
		args.update(self.primer3_global_args)
		return args

class ProbeKwargs(BaseDesignKwargs):
	n_probes: int = 100
	min_primer_probe_tm_diff: float = 6.0
	max_primer_probe_tm_diff: float = 8.0
	# Probe는 기본값이 다를 수 있으므로 재정의 혹은 상속 사용
	opt_length: int = 25
	min_length: int = 20
	max_length: int = 30
	min_gc: float = 35.0
	max_gc: float = 65.0

class BisulfiteKwargs(BaseModel):
	run: bool = False
	cpg_default: str = "methyl"

class ASPCRKwargs(BaseModel):
	cpg_default: str = "methyl"

class PCRParams(BaseModel):
	primer_kwargs: PrimerKwargs
	bisulfite_kwargs: BisulfiteKwargs
	as_pcr_kwargs: Optional[ASPCRKwargs] = None
	probe_kwargs: ProbeKwargs