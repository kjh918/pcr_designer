from __future__ import annotations

from typing import Any, Dict, List, Optional, Literal
from dataclasses import dataclass

import primer3
from Bio.Seq import reverse_complement
from Bio.SeqUtils import gc_fraction

# ==============================================================================
# 1. Genomic Region
# ==============================================================================
@dataclass
class GenomicRegion:
    chrom: str
    start: int
    end: int
    strand: str

    def to_dict(self, prefix: str = "") -> Dict[str, Any]:
        p = f"{prefix}_" if prefix else ""
        return {
            f"{p}chrom": self.chrom,
            f"{p}start": self.start,
            f"{p}end": self.end,
            f"{p}strand": self.strand,
        }

# ==============================================================================
# 2. Primer
# ==============================================================================
class Primer:
    def __init__(
        self,
        sequence: str,
        primer_type: str, # "forward", "reverse", "probe"
        
        # 위치 정보 (Optional)
        region: Optional[GenomicRegion] = None,
        relative_start: Optional[int] = None,
        relative_end: Optional[int] = None,
        
        # 실험 환경
        salt_monovalent_conc: float = 50.0,
        salt_divalent_conc: float = 1.5,
        dntp_conc: float = 0.6,
        dna_conc: float = 50.0,
        
        # 메타데이터
        penalty: float = 0.0,
    ) -> None:
        self.sequence = sequence
        self.length = len(sequence)
        self.primer_type = primer_type
        
        self.region = region
        self.relative_start = relative_start
        self.relative_end = relative_end
        
        self.salt_monovalent_conc = salt_monovalent_conc
        self.salt_divalent_conc = salt_divalent_conc
        self.dntp_conc = dntp_conc
        self.dna_conc = dna_conc
        self.penalty = penalty

        # 생성 시 자동 계산
        self._calculate_properties()

    def _calculate_properties(self) -> None:
        """예외 처리 없이 계산 수행"""
        # 1. Tm & GC
        self.tm = primer3.calc_tm(
            self.sequence,
            mv_conc=self.salt_monovalent_conc,
            dv_conc=self.salt_divalent_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
        )
        self.gc_percent = gc_fraction(self.sequence, ambiguous="ignore") * 100.0

        # 2. Secondary Structures
        args = {
            "mv_conc": self.salt_monovalent_conc,
            "dv_conc": self.salt_divalent_conc,
            "dntp_conc": self.dntp_conc,
            "dna_conc": self.dna_conc
        }
        
        hp = primer3.calc_hairpin(self.sequence, **args)
        self.hairpin_tm = hp.tm
        self.hairpin_dg = hp.dg

        hd = primer3.calc_homodimer(self.sequence, **args)
        self.homodimer_tm = hd.tm
        self.homodimer_dg = hd.dg

    # --------------------------------------------------------------------------
    # [Factory Method]
    # --------------------------------------------------------------------------
    @classmethod
    def create_from_primer3(
        cls, 
        result: Dict[str, Any], 
        rank: int, 
        role: Literal["LEFT", "RIGHT", "INTERNAL"],
        template_region: Optional[GenomicRegion] = None,
    ) -> Optional['Primer']:
        
        prefix = f"PRIMER_{role}_{rank}"
        seq = result.get(f"{prefix}_SEQUENCE")
        if not seq: return None
        
        penalty = float(result.get(f"{prefix}_PENALTY", 0.0))

        # 1. 좌표 파싱
        location_info = result.get(prefix)
        p3_index, length = location_info if location_info else (0, len(seq))

        # 2. 상대 좌표 통일
        if role == "RIGHT":
            relative_start = p3_index - length + 1
        else:
            relative_start = p3_index
        
        relative_end = relative_start + length

        # 3. 절대 좌표 계산
        region = None
        if template_region:
            if template_region.strand != "-":
                genomic_start = template_region.start + relative_start
                genomic_end = template_region.start + relative_end
            else:
                genomic_start = template_region.end - relative_end
                genomic_end = template_region.end - relative_start
            
            region = GenomicRegion(
                chrom=template_region.chrom,
                start=genomic_start,
                end=genomic_end,
                strand=template_region.strand
            )

        type_map = {"LEFT": "forward", "RIGHT": "reverse", "INTERNAL": "probe"}

        return cls(
            sequence=seq,
            primer_type=type_map.get(role, "unknown"),
            penalty=penalty,
            region=region,
            relative_start=relative_start,
            relative_end=relative_end
        )

    # --------------------------------------------------------------------------
    # Utils
    # --------------------------------------------------------------------------
    def check_gc_clamp(self) -> int:
        suffix = self.sequence[-5:]
        return suffix.count("G") + suffix.count("C")

    def count_cpg(self, template_sequence: str) -> int:
        if not template_sequence or self.relative_start is None: return 0
        n = len(template_sequence)
        
        if self.primer_type == "forward":
            start, end = max(self.relative_start, 0), min(self.relative_end + 2, n)
            return template_sequence[start:end].count("CG")
        elif self.primer_type == "reverse":
            start, end = max(self.relative_start - 1, 0), min(self.relative_end, n)
            return reverse_complement(template_sequence[start:end]).count("CG")
        return 0

    def to_dict(self) -> Dict[str, Any]:
        prefix = self.primer_type
        d = {
            f"{prefix}_sequence": self.sequence,
            f"{prefix}_length": self.length,
            f"{prefix}_tm": self.tm,
            f"{prefix}_gc": self.gc_percent,
            f"{prefix}_penalty": self.penalty,
            f"{prefix}_hairpin_tm": self.hairpin_tm,
            f"{prefix}_homodimer_tm": self.homodimer_tm,
            f"{prefix}_gc_clamp": self.check_gc_clamp()
        }
        if self.region:
            d.update(self.region.to_dict(prefix=prefix))
        return d


# ==============================================================================
# 3. Probe
# ==============================================================================
class Probe(Primer):
    def __init__(self, *args, variant_id: Optional[str] = None, allele: str = "ref", **kwargs):
        super().__init__(*args, **kwargs)
        self.variant_id = variant_id
        self.allele = allele

    @classmethod
    def create_from_primer3(cls, result, rank, template_region=None, variant_id=None, allele="ref", **kwargs):
        p = super().create_from_primer3(result, rank, "INTERNAL", template_region)
        if p:
            return cls(
                sequence=p.sequence, primer_type="probe",
                penalty=p.penalty,
                region=p.region, relative_start=p.relative_start, relative_end=p.relative_end,
                variant_id=variant_id, allele=allele
            )
        return None

    def to_dict(self) -> Dict[str, Any]:
        d = super().to_dict()
        d["probe_allele"] = self.allele
        d["probe_variant_id"] = self.variant_id
        return d


# ==============================================================================
# 4. Amplicon
# ==============================================================================
class Amplicon:
    def __init__(
        self,
        forward_primer: Primer,
        reverse_primer: Primer,
        probe: Optional[Probe] = None,
        template_sequence: str = "", 
        reference_sequence: str = "",
        pair_penalty: float = 0.0,
    ) -> None:
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.probe = probe
        self.template_sequence = template_sequence
        self.reference_sequence = reference_sequence
        self.pair_penalty = pair_penalty
        
        self.product_size = 0
        self.tm = 0.0
        self.sequence = ""
        self.region = None
        
        self._calculate_props()

    def _calculate_props(self):
        # 1. Product Size & Sequence
        if self.forward_primer.relative_start is not None and self.reverse_primer.relative_end is not None:
            self.product_size = self.reverse_primer.relative_end - self.forward_primer.relative_start
            
            if self.template_sequence:
                start = self.forward_primer.relative_start
                end = self.reverse_primer.relative_end
                if 0 <= start < end <= len(self.template_sequence):
                    self.sequence = self.template_sequence[start:end]
                    # Tm 계산 (예외처리 제거)
                    self.tm = primer3.calc_tm(self.sequence, mv_conc=50, dv_conc=1.5, dntp_conc=0.6, dna_conc=50)

        # 2. Genomic Region
        f_reg = self.forward_primer.region
        r_reg = self.reverse_primer.region
        if f_reg and r_reg:
            coords = [f_reg.start, f_reg.end, r_reg.start, r_reg.end]
            self.region = GenomicRegion(
                chrom=f_reg.chrom,
                start=min(coords),
                end=max(coords),
                strand=f_reg.strand
            )

    def to_dict(self) -> Dict[str, Any]:
        d = {
            "amplicon_sequence": self.sequence,
            "product_size": self.product_size,
            "product_tm": self.tm,
            "pair_penalty": self.pair_penalty,
        }
        if self.region: d.update(self.region.to_dict(prefix="amplicon"))
        d.update(self.forward_primer.to_dict())
        d.update(self.reverse_primer.to_dict())
        if self.probe: d.update(self.probe.to_dict())
        return d