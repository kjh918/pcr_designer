from __future__ import annotations
from typing import Any, Dict, Optional
import primer3
from Bio.SeqUtils import gc_fraction
from .primer import Primer, Probe

class Amplicon:
    """
    Forward, Reverse Primer 및 Probe 정보를 포함하는 PCR 산물 클래스.
    """
    def __init__(
        self,
        template_sequence: str,
        target_start_index: int,
        target_end_index: int,
        forward_primer: Primer,
        reverse_primer: Primer,
        probe: Optional[Probe] = None,
        product_size: Optional[int] = None
    ) -> None:
        self.template_sequence = template_sequence
        self.target_start_index = target_start_index
        self.target_end_index = target_end_index
        
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.probe = probe
        
        # 1. Amplicon 범위 및 서열 추출
        # Fwd의 시작점(5')부터 Rev의 끝점(3'의 템플릿 상 위치)까지
        # 주의: Primer.binding_end_index가 Python slicing end(exclusive) 기준이라고 가정
        if (self.forward_primer.binding_start_index is not None and 
            self.reverse_primer.binding_end_index is not None):
            
            self.start_index = self.forward_primer.binding_start_index
            self.end_index = self.reverse_primer.binding_end_index
            
            self.sequence = self.template_sequence[self.start_index : self.end_index]
            self.product_size = len(self.sequence)
        else:
            self.sequence = ""
            self.product_size = 0
            
        # 외부에서 주입된 사이즈와 계산된 사이즈가 다르면 로깅/워닝 가능
        if product_size and product_size != self.product_size:
            # print(f"Warning: Calculated size ({self.product_size}) differs from input ({product_size})")
            pass

        # 2. 물성 계산
        self.tm = 0.0
        self.gc_percent = 0.0
        if self.sequence:
            self._calc_properties()

    def _calc_properties(self):
        # Primer와 동일한 농도 조건 사용
        self.tm = primer3.calc_tm(
            self.sequence,
            mv_conc=self.forward_primer.salt_monovalent_conc,
            dv_conc=self.forward_primer.salt_divalent_conc,
            dntp_conc=self.forward_primer.dntp_conc,
            dna_conc=self.forward_primer.dna_conc,
        )
        self.gc_percent = gc_fraction(self.sequence, ambiguous="ignore") * 100.0

    def to_dict(self) -> Dict[str, Any]:
        """
        Amplicon 및 포함된 Primer들의 정보를 Flat Dictionary로 변환 (리포트용)
        """
        data = {
            "amplicon_sequence": self.sequence,
            "product_size": self.product_size,
            "product_tm": self.tm,
            "product_gc": self.gc_percent,
            "target_start": self.target_start_index,
            "target_end": self.target_end_index
        }
        
        # Fwd, Rev Primer 정보 병합
        data.update(self.forward_primer.to_dict())
        data.update(self.reverse_primer.to_dict())
        
        # Probe 정보 병합 (있을 경우)
        if self.probe:
            data.update(self.probe.to_dict())
            
        return data