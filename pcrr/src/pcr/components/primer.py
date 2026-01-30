from __future__ import annotations
from .targets import TargetType, Variant, CpG 
from typing import Any, Dict, List, Optional, Literal
from pydantic import BaseModel, Field

import primer3
from Bio.Seq import reverse_complement
from Bio.SeqUtils import gc_fraction
from ..utils import get_start_end_index

Allele = Literal["ref", "alt"]

class Primer(BaseModel):
    # --------------------------------------------------------
    # 1. 필드 선언 (Pydantic이 인식하도록 상단에 배치)
    # --------------------------------------------------------
    template_sequence: str
    sequence: str
    strand: str
    primer_type: str
    target_start_index: int
    target_end_index: int
        
    # 선택적 필드
    reference_template_sequence: Optional[str] = None
    chrom: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    binding_start_index: Optional[int] = None
    binding_end_index: Optional[int] = None
        
    # Tm
    tm: Optional[float] = Field(None, description="Melting Temperature of the primer")
        
    # 설정값
    salt_monovalent_conc: float = 50.0
    salt_divalent_conc: float = 1.5
    dntp_conc: float = 0.6
    dna_conc: float = 50.0

    # 계산된 필드 (미리 선언 필수)
    length: int = 0
    gc_percent: Optional[float] = None
    
    hairpin: bool = False
    hairpin_tm: Optional[float] = None
    hairpin_dg: Optional[float] = None
    
    homodimer: bool = False
    homodimer_tm: Optional[float] = None
    homodimer_dg: Optional[float] = None

    # --------------------------------------------------------
    # 2. __init__ 메서드 (로직을 함수 안으로 넣음)
    # --------------------------------------------------------
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
        tm: Optional[float] = None,
        salt_monovalent_conc: float = 50.0,
        salt_divalent_conc: float = 1.5,
        dntp_conc: float = 0.6,
        dna_conc: float = 50.0,
        **kwargs # Pydantic 호환용
    ) -> None:
        # [중요] Pydantic의 초기화 로직을 먼저 실행
        super().__init__(
            template_sequence=template_sequence,
            sequence=sequence,
            strand=strand,
            primer_type=primer_type,
            target_start_index=target_start_index,
            target_end_index=target_end_index,
            reference_template_sequence=reference_template_sequence,
            chrom=chrom,
            start=start,
            end=end,
            binding_start_index=binding_start_index,
            binding_end_index=binding_end_index,
            tm=tm,
            salt_monovalent_conc=salt_monovalent_conc,
            salt_divalent_conc=salt_divalent_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
            **kwargs
        )

        # [중요] 계산 로직은 super().__init__ 이후에 실행
        self.length = len(sequence)
        
        # 좌표 설정 로직 (값이 안 들어왔을 때만 계산)
        if self.binding_start_index is None or self.binding_end_index is None:
             self.binding_start_index, self.binding_end_index = get_start_end_index(self.template_sequence, self.sequence)
        
        # 물성 계산 메서드 호출
        self._calc_basic_properties()
        self._calc_hairpin()
        self._calc_homodimer()

    # --------------------------------------------------------
    # 3. 내부 계산 메서드
    # --------------------------------------------------------
    def _calc_basic_properties(self) -> None:
        if self.tm is None:
            try:
                self.tm = primer3.calc_tm(
                    self.sequence,
                    mv_conc=self.salt_monovalent_conc,
                    dv_conc=self.salt_divalent_conc,
                    dntp_conc=self.dntp_conc,
                    dna_conc=self.dna_conc,
                )
            except Exception:
                self.tm = 0.0
        self.gc_percent = gc_fraction(self.sequence, ambiguous="ignore") * 100.0

    def _calc_hairpin(self) -> None:
        try:
            result = primer3.calc_hairpin(
                self.sequence,
                mv_conc=self.salt_monovalent_conc,
                dv_conc=self.salt_divalent_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
            )
            self.hairpin = result.structure_found
            self.hairpin_tm = result.tm
            self.hairpin_dg = result.dg / 1000.0 if result.dg else 0.0
        except Exception:
            pass

    def _calc_homodimer(self) -> None:
        try:
            result = primer3.calc_homodimer(
                self.sequence,
                mv_conc=self.salt_monovalent_conc,
                dv_conc=self.salt_divalent_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
            )
            self.homodimer = result.structure_found
            self.homodimer_tm = result.tm
            self.homodimer_dg = result.dg / 1000.0 if result.dg else 0.0
        except Exception:
            pass

    # -------------------------
    # checks
    # -------------------------
    def check_three_prime_is(self, sequence: str, *, use_reference: bool = True) -> bool:
        ref = self.reference_template_sequence if use_reference else self.template_sequence
        seq_len = len(sequence)

        if self.strand == "forward":
            # end_index가 None일 경우 안전장치 필요
            if self.binding_end_index is None: return False
            end = self.binding_end_index 
            start = max(end - seq_len, 0)
            three_prime_seq = ref[start:end]
        elif self.strand == "reverse":
            if self.binding_start_index is None: return False
            start = self.binding_start_index
            end = min(self.binding_start_index + seq_len, len(ref))
            three_prime_seq = reverse_complement(ref[start:end])
        else:
            raise ValueError(f"Unknown strand type: {self.strand}")

        return three_prime_seq == sequence
        
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
        # Pydantic v2 호환
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
    # Probe 전용 필드 추가
    target: Optional[TargetType] = None

    def __init__(self, *args, target: Optional[TargetType] = None, **kwargs) -> None:
        # primer_type 강제 지정
        kwargs["primer_type"] = "probe"
        super().__init__(*args, **kwargs)
        self.target = target

    def covers_target(self) -> bool:
        """Probe가 Target의 전체 구간을 포함하는지 확인"""
        if not self.target:
            return False
        
        if self.binding_start_index is None or self.binding_end_index is None:
            return False

        return (self.binding_start_index <= self.target.index_start) and \
               (self.target.index_end <= self.binding_end_index)

    def validate_specificity(self) -> bool:
        if not self.covers_target():
            return False

        expected_seq = self.target.get_expected_sequence().upper()

        rel_start = self.target.index_start - self.binding_start_index
        rel_end = self.target.index_end - self.binding_start_index
        
        if rel_start < 0 or rel_end > len(self.sequence):
            return False
            
        actual_seq = self.sequence[rel_start:rel_end].upper()

        return actual_seq == expected_seq

    def to_dict(self, ignore_attributes: Optional[List[str]] = None) -> Dict[str, Any]:
        d = super().to_dict(ignore_attributes)
        
        if self.target:
            d["target_type"] = self.target.__class__.__name__
            d["target_chrom"] = getattr(self.target, 'chrom', None)
            
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