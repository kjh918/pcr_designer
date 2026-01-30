from __future__ import annotations
from .targets import TargetType, Variant, CpG 
from typing import Any, Dict, List, Optional, Literal
from pydantic import BaseModel, Field

import primer3
from Bio.Seq import reverse_complement
from Bio.SeqUtils import gc_fraction
from ..utils import get_start_end_index

# ----------------------------------------------------------------------
# 3. Amplicon Class (Standard Python Class)
# ----------------------------------------------------------------------
class Amplicon:
    """
    Forward, Reverse Primer 및 Probe 정보를 포함하는 PCR 산물 클래스.
    """
    def __init__(
        self,
        template_sequence: str,
        target_start_index: int,
        target_end_index: int,
        forward_primer: Optional[Primer] = None, # Probe Only 단계 때문에 Optional 처리
        reverse_primer: Optional[Primer] = None,
        probe: Optional[Probe] = None,
        product_size: Optional[int] = None,
        tm: Optional[float] = None,
        reference_template_sequence: Optional[str] = None,
        # QC 결과 저장을 위한 필드 추가
        qc_status: Optional[Any] = None,
        is_qc_pass: bool = False,
        off_target_count: int = 0
    ) -> None:
        # 1. 입력값 저장
        self.template_sequence = template_sequence
        self.target_start_index = target_start_index
        self.target_end_index = target_end_index
        
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.probe = probe
        self.reference_template_sequence = reference_template_sequence
        
        # QC 상태 저장
        self.qc_status = qc_status
        self.is_qc_pass = is_qc_pass
        self.off_target_count = off_target_count

        # 2. Amplicon 범위 및 서열 추출 로직
        self.start_index: Optional[int] = None
        self.end_index: Optional[int] = None
        self.sequence: str = ""
        self.product_size: int = 0
        
        # Fwd와 Rev가 모두 존재하고, 바인딩 좌표가 있을 때만 계산
        if (self.forward_primer and self.forward_primer.binding_start_index is not None and 
            self.reverse_primer and self.reverse_primer.binding_end_index is not None):
            
            self.start_index = self.forward_primer.binding_start_index
            self.end_index = self.reverse_primer.binding_end_index
            
            # 슬라이싱 (Python은 end_index가 exclusive이므로 그대로 사용)
            self.sequence = self.template_sequence[self.start_index : self.end_index]
            self.product_size = len(self.sequence)
        
        # 외부(Primer3)에서 계산된 사이즈가 들어왔다면 우선순위로 덮어쓰거나 검증 가능
        if product_size is not None:
             # 만약 계산된 것과 다르면 Primer3 값을 신뢰하거나 로깅
             self.product_size = product_size

        # 3. 물성 계산 (Tm, GC)
        self.tm = tm # 외부 입력값 우선
        self.gc_percent: float = 0.0
        
        if self.sequence:
            self._calc_properties()

    def _calc_properties(self) -> None:
        """Tm 및 GC 함량 계산"""
        # Tm이 없을 경우에만 계산
        if self.tm is None and self.forward_primer:
            try:
                # Primer와 동일한 Salt 조건 사용
                self.tm = primer3.calc_tm(
                    self.sequence,
                    mv_conc=self.forward_primer.salt_monovalent_conc,
                    dv_conc=self.forward_primer.salt_divalent_conc,
                    dntp_conc=self.forward_primer.dntp_conc,
                    dna_conc=self.forward_primer.dna_conc,
                )
            except Exception:
                self.tm = 0.0
        
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
            "target_end": self.target_end_index,
            "qc_pass": self.is_qc_pass,
            "off_targets": self.off_target_count
        }
        
        # Fwd, Rev Primer 정보 병합 (있을 경우만)
        if self.forward_primer:
            data.update(self.forward_primer.to_dict())
        if self.reverse_primer:
            data.update(self.reverse_primer.to_dict())
        
        # Probe 정보 병합 (있을 경우)
        if self.probe:
            data.update(self.probe.to_dict())
            
        return data