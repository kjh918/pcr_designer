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
    end: int    # 0-based end (exclusive)
    strand: str # "+" or "-"
    sequence: str
	base: int = 0
	
    def __len__(self) -> int:
        return self.end - self.start

    def to_dict(self, prefix: str = "") -> Dict[str, Any]:
        """Flatten하여 반환 (예: forward_primer_chrom)"""
        p = f"{prefix}_" if prefix else ""
        return {
            f"{p}chrom": self.chrom,
            f"{p}start": self.start,
            f"{p}end": self.end,
            f"{p}strand": self.strand,
            f"{p}base": self.base,
        }
		
    def to_dict(self, prefix: str = "") -> Dict[str, Any]:
        """Flatten하여 반환 (예: forward_primer_chrom)"""
        p = f"{prefix}_" if prefix else ""
        return {
            f"{p}chrom": self.chrom,
            f"{p}start": self.start,
            f"{p}end": self.end,
            f"{p}strand": self.strand,
        }