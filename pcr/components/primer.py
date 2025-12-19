from __future__ import annotations

from typing import Any, Dict, List, Optional, Literal

import primer3
from Bio.Seq import reverse_complement
from Bio.SeqUtils import gc_fraction

from pcr.utils import get_start_end_index
from pcr.components.variant import Variant

Allele = Literal["ref", "alt"]

class Primer:
    # ... (기존 필드들 동일)

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
        # ✅ 추가: mismatch primer는 템플릿에서 검색이 안되므로 바인딩 좌표를 직접 주입
        binding_start_index: Optional[int] = None,
        binding_end_index: Optional[int] = None,
        # ... (salt/dntp/dna_conc 동일)
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

        print(self.primer_type, self.template_sequence, self.sequence)

        # ✅ mismatch primer 대응
        if binding_start_index is not None and binding_end_index is not None:
            self.start_index = binding_start_index
            self.end_index = binding_end_index
        else:
            print(self.template_sequence, self.sequence)
            self.start_index, self.end_index = get_start_end_index(self.template_sequence, self.sequence)

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

    # -------------------------
    # thermo
    # -------------------------
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
    # checks (reference/template 선택)
    # -------------------------
    def check_three_prime_is(self, sequence: str, *, use_reference: bool = True) -> bool:
        ref = self.reference_template_sequence if use_reference else self.template_sequence
        seq_len = len(sequence)

        if self.strand == "forward":
            end = self.end_index + 1
            start = max(end - seq_len, 0)
            three_prime_seq = ref[start:end]
        elif self.strand == "reverse":
            start = self.start_index
            end = min(self.start_index + seq_len, len(ref))
            three_prime_seq = reverse_complement(ref[start:end])
        else:
            raise ValueError(f"Unknown strand type: {self.strand}")

        return three_prime_seq == sequence

    def count_cpg(self, *, use_reference: bool = True) -> int:
        ref = self.reference_template_sequence if use_reference else self.template_sequence
        n = len(ref)

        if self.strand == "forward":
            start = max(self.start_index, 0)
            end = min(self.end_index + 2, n)
            return ref[start:end].count("CG")

        if self.strand == "reverse":
            start = max(self.start_index - 1, 0)
            end = min(self.end_index + 1, n)
            window = reverse_complement(ref[start:end])
            return window.count("CG")

        raise ValueError(f"Unknown strand type: {self.strand}")

    def count_non_cpg_cytosine(self, *, use_reference: bool = True) -> int:
        ref = self.reference_template_sequence if use_reference else self.template_sequence

        if self.strand == "forward":
            window = ref[self.start_index : self.end_index + 1]
        elif self.strand == "reverse":
            window = reverse_complement(ref[self.start_index : self.end_index + 1])
        else:
            raise ValueError(f"Unknown strand type: {self.strand}")

        total_c = window.count("C")
        return total_c - self.count_cpg(use_reference=use_reference)

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
    Single-variant targeting probe.
    reference_template_sequence는 ref allele 기준(좌표 고정),
    template_sequence는 ref/alt + conversion 반영(allele-specific 검증은 template 기준).
    """
    def __init__(self, *args, variant: Variant, **kwargs) -> None:
        super().__init__(*args, primer_type="probe", **kwargs)
        self.variant = variant

    def covers_variant(self) -> bool:
        return self.start_index <= self.variant.index <= self.end_index

    def allele_base_at_variant(self, *, use_reference: bool = False) -> str:
        seq = self.reference_template_sequence if use_reference else self.template_sequence
        return seq[self.variant.index]

    def validate_variant_specific(self) -> bool:
        if not self.covers_variant():
            return False
        expected = self.variant.ref.upper() if self.allele == "ref" else self.variant.alt.upper()
        actual = self.allele_base_at_variant(use_reference=False).upper()
        return actual == expected
