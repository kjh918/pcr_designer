# designers/as_pcr.py
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import primer3

from pcr.components.primer import Primer
from pcr.components.as_assay import AsPcrAssay


DNA_BASES = ("A", "C", "G", "T")


def _mismatch_bases(template_base: str) -> List[str]:
    b = template_base.upper()
    return [x for x in DNA_BASES if x != b]


def _build_as_forward_variants(
    *,
    reference_template_sequence: str,
    template_sequence_for_calc: str,
    target_start_index: int,
    ref_allele: str,
    alt_allele: str,
    min_len: int,
    max_len: int,
    mismatch_offset_from_3p: int = 3,  # 3'에서 3번째 = offset 3
) -> List[Tuple[str, str, int, int]]:
    """
    반환: (wt_seq, alt_seq, binding_start, binding_end) 후보들
    - binding_end는 항상 target_start_index
    - mismatch는 end 기준으로 (end-(offset-1)) 위치에 넣음
      예: offset=3 => end-2 위치가 mismatch
    """
    ref = reference_template_sequence
    n = len(ref)
    end = target_start_index
    if not (0 <= end < n):
        return []

    # 템플릿의 변이 위치 기준 염기(레퍼런스)
    template_3p_base = ref[end].upper()
    if template_3p_base != ref_allele.upper():
        # reference_template_sequence가 ref allele 기준이라는 전제 위반
        # (너가 말한 규칙과 다르면 여기서 early stop)
        return []

    out: List[Tuple[str, str, int, int]] = []

    mismatch_pos = end - (mismatch_offset_from_3p - 1)  # offset=3 => end-2
    if mismatch_pos < 0:
        return []

    template_mismatch_base = ref[mismatch_pos].upper()
    mismatch_choices = _mismatch_bases(template_mismatch_base)

    for L in range(min_len, max_len + 1):
        start = end - L + 1
        if start < 0:
            continue

        # 기본 프라이머 서열(레퍼런스 템플릿에서 뽑음)
        core = list(ref[start : end + 1].upper())
        if len(core) != L:
            continue

        # 3' 말단 고정 (WT=ref, ALT=alt)
        wt = core[:]
        alt = core[:]
        wt[-1] = ref_allele.upper()
        alt[-1] = alt_allele.upper()

        # -3 mismatch 적용 (모든 가능한 mismatch base로 후보 확장)
        idx_in_primer = mismatch_pos - start  # primer 내 index
        if not (0 <= idx_in_primer < L):
            continue

        for mm in mismatch_choices:
            wt2 = wt[:]
            alt2 = alt[:]
            wt2[idx_in_primer] = mm
            alt2[idx_in_primer] = mm

            out.append(("".join(wt2), "".join(alt2), start, end))

    return out


class AsPcrDesigner:
    """
    AS-PCR: forward 2개(wt/alt) + reverse 1개
    - forward 3' end == target_start_index
    - WT 3' == ref_allele
    - ALT 3' == alt_allele
    - 3'에서 3번째(-3) 염기 mismatch
    """

    def __init__(
        self,
        *,
        template_sequence: str,
        reference_template_sequence: str,
        target_start_index: int,
        target_end_index: int,
        ref_allele: str,
        alt_allele: str,
        min_amplicon_length: int,
        max_amplicon_length: int,
        n_reverse: int = 50,
        forward_min_len: int = 18,
        forward_max_len: int = 28,
        primer3_global_args: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.template_sequence = template_sequence
        self.reference_template_sequence = reference_template_sequence
        self.target_start_index = target_start_index
        self.target_end_index = target_end_index

        self.ref_allele = ref_allele
        self.alt_allele = alt_allele

        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length

        self.n_reverse = n_reverse
        self.forward_min_len = forward_min_len
        self.forward_max_len = forward_max_len

        self.primer3_global_args = primer3_global_args or {}

        self.assays: List[AsPcrAssay] = []

    def reset(self) -> None:
        self.assays = []

    def _design_reverse_primers(self) -> List[Tuple[str, int, int]]:
        """
        primer3로 reverse primer 후보를 뽑는다.
        return: (rev_seq, rev_start, rev_len)
        """
        seq_args = {
            "SEQUENCE_ID": "AS_PCR",
            "SEQUENCE_TEMPLATE": self.template_sequence,
            # target은 그대로 두되, reverse primer만 뽑는 용도
            "SEQUENCE_TARGET": [self.target_start_index, self.target_end_index - self.target_start_index + 1],
        }
        global_args = {
            "PRIMER_TASK": "generic",
            "PRIMER_NUM_RETURN": self.n_reverse,
            "PRIMER_PICK_LEFT_PRIMER": 0,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 0,
            "PRIMER_PRODUCT_SIZE_RANGE": [self.min_amplicon_length, self.max_amplicon_length],
        }
        global_args.update(self.primer3_global_args)

        res = primer3.bindings.designPrimers(seq_args=seq_args, global_args=global_args)

        out: List[Tuple[str, int, int]] = []
        n = res.get("PRIMER_RIGHT_NUM_RETURNED", 0) or 0
        for i in range(n):
            seq = res.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
            loc = res.get(f"PRIMER_RIGHT_{i}")  # [start, len]
            if not seq or not loc:
                continue
            start, ln = int(loc[0]), int(loc[1])
            out.append((seq, start, ln))
        return out

    def design(self) -> List[AsPcrAssay]:
        self.assays = []

        reverse_candidates = self._design_reverse_primers()

        # forward 후보(서열/좌표) 생성: reference 기준으로 뽑고, WT/ALT 생성
        fw_variants = _build_as_forward_variants(
            reference_template_sequence=self.reference_template_sequence,
            template_sequence_for_calc=self.template_sequence,
            target_start_index=self.target_start_index,
            ref_allele=self.ref_allele,
            alt_allele=self.alt_allele,
            min_len=self.forward_min_len,
            max_len=self.forward_max_len,
            mismatch_offset_from_3p=3,
        )

        # reverse 후보와 조합 -> amplicon size 필터
        for (rev_seq, rev_start, rev_len) in reverse_candidates:
            # primer3 right primer의 start는 "leftmost index" 기준(라이브러리 표현)
            rev_binding_start = rev_start
            rev_binding_end = rev_start + rev_len - 1

            # forward end는 target_start_index
            for (wt_seq, alt_seq, fw_start, fw_end) in fw_variants:
                # amplicon: fw_start ~ rev_binding_end (단, reverse가 forward보다 뒤에 있어야 함)
                if rev_binding_end <= fw_end:
                    continue

                amp_len = rev_binding_end - fw_start + 1
                if not (self.min_amplicon_length <= amp_len <= self.max_amplicon_length):
                    continue

                # Primer 객체 생성 (mismatch 때문에 binding index 주입!)
                wt_fw = Primer(
                    template_sequence=self.template_sequence,
                    reference_template_sequence=self.reference_template_sequence,
                    sequence=wt_seq,
                    strand="forward",
                    primer_type="wt_forward",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=fw_start,
                    binding_end_index=fw_end,
                )
                alt_fw = Primer(
                    template_sequence=self.template_sequence,
                    reference_template_sequence=self.reference_template_sequence,
                    sequence=alt_seq,
                    strand="forward",
                    primer_type="alt_forward",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=fw_start,
                    binding_end_index=fw_end,
                )
                rev = Primer(
                    template_sequence=self.template_sequence,
                    reference_template_sequence=self.reference_template_sequence,
                    sequence=rev_seq,
                    strand="reverse",
                    primer_type="reverse",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    # reverse는 템플릿에 그대로 존재하므로 binding index 없어도 됨
                )

                amplicon_seq = self.template_sequence[fw_start : rev_binding_end + 1]

                self.assays.append(
                    AsPcrAssay(
                        wt_forward=wt_fw,
                        alt_forward=alt_fw,
                        reverse=rev,
                        amplicon_start_index=fw_start,
                        amplicon_end_index=rev_binding_end,
                        amplicon_sequence=amplicon_seq,
                    )
                )

        return self.assays
