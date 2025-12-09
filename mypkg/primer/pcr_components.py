from __future__ import annotations

from typing import Dict, Any, List, Optional, Tuple

import primer3
from Bio.Seq import reverse_complement
from Bio.SeqUtils import gc_fraction


# ----------------------------------------------------------------------
# 유틸 함수
# ----------------------------------------------------------------------
def get_start_end_index(template_sequence: str, sequence: str) -> Tuple[int, int]:
    """
    template_sequence에서 주어진 sequence(또는 그 reverse complement)의
    시작/끝 인덱스를 반환한다. (0-based, end inclusive)

    찾지 못하면 ValueError 발생.
    """
    try:
        start_index = template_sequence.index(sequence)
    except ValueError:
        rc = reverse_complement(sequence)
        try:
            start_index = template_sequence.index(rc)
        except ValueError:
            raise ValueError(f"{sequence} or its reverse complement not found in template.")

    end_index = start_index + len(sequence) - 1
    return start_index, end_index


# ----------------------------------------------------------------------
# Primer
# ----------------------------------------------------------------------
class Primer:
    template_sequence: str
    sequence: str
    strand: str  # 'forward' or 'reverse'
    primer_type: str  # 'forward' / 'reverse' / 'probe'
    tm: float
    gc_percent: float

    DEFAULT_SALT_MONOVALENT: float = 50.0
    DEFAULT_SALT_DIVALENT: float = 1.5
    DEFAULT_DNTP_CONC: float = 0.6
    DEFAULT_DNA_CONC: float = 50.0

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
        salt_monovalent_conc: float = DEFAULT_SALT_MONOVALENT,
        salt_divalent_conc: float = DEFAULT_SALT_DIVALENT,
        dntp_conc: float = DEFAULT_DNTP_CONC,
        dna_conc: float = DEFAULT_DNA_CONC,
    ) -> None:
        # 기본 서열/위치 정보
        self.template_sequence: str = template_sequence
        self.reference_template_sequence: str = reference_template_sequence or template_sequence
        self.sequence: str = sequence
        self.strand: str = strand
        self.primer_type: str = primer_type
        self.target_start_index: int = target_start_index
        self.target_end_index: int = target_end_index

        self.length: int = len(sequence)
        self.start_index, self.end_index = get_start_end_index(
            self.template_sequence, self.sequence
        )

        # 위치 정보(게놈 상)
        self.chrom: Optional[str] = chrom
        self.start: Optional[int] = start
        self.end: Optional[int] = end

        # PCR 반응 조건
        self.salt_monovalent_conc: float = salt_monovalent_conc
        self.salt_divalent_conc: float = salt_divalent_conc
        self.dntp_conc: float = dntp_conc
        self.dna_conc: float = dna_conc

        # Thermodynamics / 특성 계산
        self._calc_basic_properties()
        self._calc_hairpin()
        self._calc_homodimer()

    # ------------------------------------------------------------------
    # 내부 계산 메서드
    # ------------------------------------------------------------------
    def _calc_basic_properties(self) -> None:
        """Tm, GC% 등 기본 특성을 계산."""
        self.tm = primer3.calc_tm(
            self.sequence,
            mv_conc=self.salt_monovalent_conc,
            dv_conc=self.salt_divalent_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
        )
        self.gc_percent = gc_fraction(self.sequence, ambiguous="ignore") * 100.0

    def _calc_hairpin(self) -> None:
        """Hairpin 구조 가능성 계산."""
        result = primer3.calc_hairpin(
            self.sequence,
            mv_conc=self.salt_monovalent_conc,
            dv_conc=self.salt_divalent_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
        )
        self.hairpin: bool = result.structure_found
        self.hairpin_tm: float = result.tm
        # primer3 ThermoResult는 dg/dh/ds를 1000배한 값으로 주기 때문에 다시 환산
        self.hairpin_dg: float = result.dg / 1000.0
        self.hairpin_dh: float = result.dh / 1000.0
        self.hairpin_ds: float = result.ds / 1000.0

    def _calc_homodimer(self) -> None:
        """Homodimer 구조 가능성 계산."""
        result = primer3.calc_homodimer(
            self.sequence,
            mv_conc=self.salt_monovalent_conc,
            dv_conc=self.salt_divalent_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
        )
        self.homodimer: bool = result.structure_found
        self.homodimer_tm: float = result.tm
        self.homodimer_dg: float = result.dg / 1000.0
        self.homodimer_dh: float = result.dh / 1000.0
        self.homodimer_ds: float = result.ds / 1000.0

    # ------------------------------------------------------------------
    # 특성 체크 메서드
    # ------------------------------------------------------------------
    def check_three_prime_is(self, sequence: str) -> bool:
        """
        3' 말단 부분이 주어진 sequence와 일치하는지 확인.
        forward: 템플릿에서 프라이머 끝 쪽
        reverse: reverse_complement 템플릿 상에서 프라이머 끝 쪽
        """
        ref = self.reference_template_sequence
        seq_len = len(sequence)

        if self.strand == "forward":
            end = self.end_index + 1
            start = max(end - seq_len, 0)
            three_prime_seq = ref[start:end]
        elif self.strand == "reverse":
            # reverse primer는 reverse complement 방향으로 3' 말단 확인
            start = self.start_index
            end = min(self.start_index + seq_len, len(ref))
            three_prime_seq = reverse_complement(ref[start:end])
        else:
            raise ValueError(f"Unknown strand type: {self.strand}")

        return three_prime_seq == sequence

    def check_cpg_count(self, min_cpg_count: int, max_cpg_count: int) -> bool:
        cpg_count = self.count_cpg()
        return min_cpg_count <= cpg_count <= max_cpg_count

    def check_non_cpg_cytosine_count(
        self, min_non_cpg_cytosine_count: int, max_non_cpg_cytosine_count: int) -> bool:
        non_cpg_c = self.count_non_cpg_cytosine()
        return min_non_cpg_cytosine_count <= non_cpg_c <= max_non_cpg_cytosine_count

    def count_cpg(self) -> int:
        """
        프라이머 binding region에서 CpG count.
        forward: [start_index : end_index+2]
        reverse: reverse_complement([start_index-1 : end_index+1])
        (인덱스가 범위를 벗어나면 자동으로 clamp)
        """
        ref = self.reference_template_sequence
        n = len(ref)

        if self.strand == "forward":
            start = max(self.start_index, 0)
            end = min(self.end_index + 2, n)
            window = ref[start:end]
            return window.count("CG")

        if self.strand == "reverse":
            start = max(self.start_index - 1, 0)
            end = min(self.end_index + 1, n)
            window = reverse_complement(ref[start:end])
            return window.count("CG")

        raise ValueError(f"Unknown strand type: {self.strand}")

    def count_non_cpg_cytosine(self) -> int:
        """
        프라이머 영역 내 non-CpG cytosine 개수.
        전체 C 개수 - CpG 개수
        """
        ref = self.reference_template_sequence
        if self.strand == "forward":
            window = ref[self.start_index : self.end_index + 1]
        elif self.strand == "reverse":
            window = reverse_complement(ref[self.start_index : self.end_index + 1])
        else:
            raise ValueError(f"Unknown strand type: {self.strand}")

        total_c = window.count("C")
        return total_c - self.count_cpg()

    def to_dict(
        self,
        ignore_attributes: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Primer 정보를 dict로 변환.
        key 이름은 "{primer_type}_{attribute}" 형식으로 붙인다.
        """
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

        primer_dict: Dict[str, Any] = {}
        for key, value in self.__dict__.items():
            if key in ignore_attributes:
                continue
            primer_dict[f"{self.primer_type}_{key}"] = value

        return primer_dict


# ----------------------------------------------------------------------
# Amplicon
# ----------------------------------------------------------------------
class Amplicon:
    reference_template_sequence: str
    template_sequence: str
    target_start_index: int
    target_end_index: int
    chrom: Optional[str]
    start: Optional[int]
    end: Optional[int]
    forward_primer: Optional[Primer]
    reverse_primer: Optional[Primer]
    probe: Optional[Primer]
    amplicon_sequence: Optional[str]

    def __init__(
        self,
        template_sequence: str,
        target_start_index: int,
        target_end_index: int,
        reference_template_sequence: Optional[str] = None,
        chrom: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        forward_primer: Optional[Primer] = None,
        reverse_primer: Optional[Primer] = None,
        probe: Optional[Primer] = None,
    ) -> None:
        self.template_sequence: str = template_sequence
        self.reference_template_sequence: str = reference_template_sequence or template_sequence

        self.target_start_index: int = target_start_index
        self.target_end_index: int = target_end_index

        self.chrom: Optional[str] = chrom
        self.start: Optional[int] = start
        self.end: Optional[int] = end

        self.forward_primer: Optional[Primer] = forward_primer
        self.reverse_primer: Optional[Primer] = reverse_primer
        self.probe: Optional[Primer] = probe

        if self.forward_primer is not None:
            self.forward_start_index, self.forward_end_index = get_start_end_index(
                self.template_sequence, self.forward_primer.sequence
            )

        if self.reverse_primer is not None:
            self.reverse_start_index, self.reverse_end_index = get_start_end_index(
                self.template_sequence, self.reverse_primer.sequence
            )

        if self.probe is not None:
            self.probe_start_index, self.probe_end_index = get_start_end_index(
                self.template_sequence, self.probe.sequence
            )

        self.amplicon_sequence: Optional[str] = self._calc_amplicon_sequence()

    def _calc_amplicon_sequence(self) -> Optional[str]:
        """forward / reverse 프라이머가 둘 다 있는 경우만 amplicon 시퀀스를 계산."""
        if self.forward_primer is None or self.reverse_primer is None:
            return None
        return self.template_sequence[self.forward_start_index : self.reverse_end_index + 1]

    def to_dict(self) -> Dict[str, Any]:
        amplicon_dict: Dict[str, Any] = {
            "reference_template_sequence": self.reference_template_sequence,
            "template_sequence": self.template_sequence,
            "target_start_index": self.target_start_index,
            "target_end_index": self.target_end_index,
        }

        if self.amplicon_sequence is not None:
            amplicon_dict["amplicon_sequence"] = self.amplicon_sequence
            amplicon_dict["amplicon_length"] = len(self.amplicon_sequence)
        else:
            amplicon_dict["amplicon_sequence"] = None
            amplicon_dict["amplicon_length"] = None

        if self.forward_primer is not None:
            amplicon_dict.update(self.forward_primer.to_dict())
        if self.reverse_primer is not None:
            amplicon_dict.update(self.reverse_primer.to_dict())
        if self.probe is not None:
            amplicon_dict.update(self.probe.to_dict())

        return amplicon_dict
