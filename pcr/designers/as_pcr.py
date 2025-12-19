# pcr/designers/as_pcr.py
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import primer3

from pcr.designers.base import BasePrimerDesigner
from pcr.components import Amplicon
from pcr.components.primer import Primer


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------


def _replace_3prime_base(primer_5to3: str, base: str, strand: str = "forward") -> str:
    """
    primer3가 준 primer 서열(항상 5'->3')에서 3' 말단 염기를 강제 치환.

    - forward primer의 3' end: 문자열 마지막(s[-1])
    - reverse primer의 3' end: 문자열 첫 글자(s[0])  (5'->3' 표기에서 왼쪽이 3' end에 해당)
      -> 따라서 reverse는 s[0]을 바꾼다.
    """
    change_info_dict = {
        'A':'T',
        'T':'A',
        'G':'C',
        'C':'G',
    }
    s = primer_5to3.upper()
    if not s:
        return s
    b = base.upper()

    if strand == "forward":
        return s[:-1] + b
    elif strand == "reverse":
        return s[:-1] + change_info_dict[b]
    else:
        raise ValueError(f"Unknown strand: {strand}")


def _mk_primer(
    *,
    template_sequence: str,
    reference_template_sequence: str,
    sequence: str,
    strand: str,
    primer_type: str,
    target_start_index: int,
    target_end_index: int,
    binding_start_index: Optional[int] = None,
    binding_end_index: Optional[int] = None,
) -> Primer:
    """
    프로젝트 Primer.__init__ 시그니처에 맞춘 factory.
    mismatch/치환 등으로 template에서 exact match가 안 나도
    binding_start_index/binding_end_index 주입으로 좌표를 고정한다.
    """
    return Primer(
        template_sequence=template_sequence,
        reference_template_sequence=reference_template_sequence,
        sequence=sequence,
        strand=strand,
        primer_type=primer_type,
        target_start_index=target_start_index,
        target_end_index=target_end_index,
        binding_start_index=binding_start_index,
        binding_end_index=binding_end_index,
    )


def _pos_to_binding_strict(pos: Any) -> Tuple[int, int]:
    """
    primer3의 PRIMER_LEFT_i / PRIMER_RIGHT_i = [start, length] 를
    (binding_start, binding_end)로 변환. 실패 시 예외.
    """
    if not pos or len(pos) != 2:
        raise ValueError(f"Invalid primer3 position: {pos}")
    start, length = pos
    start = int(start)
    length = int(length)
    if length <= 0:
        raise ValueError(f"Invalid primer length: {length}")
    return start, start + length - 1


# ----------------------------------------------------------------------
# AsPcrDesigner
# ----------------------------------------------------------------------
class AsPcrDesigner(BasePrimerDesigner):
    """
    AS-PCR 디자이너:
    - 입력 template_sequence: 보통 alt_template_sequence (변이 반영 서열)
    - 입력 reference_template_sequence: 보통 ref_template_sequence (ref 반영 서열)

    생성 Amplicon:
      - ref amplicon: ref_template_sequence 기반으로 Primer/Amplicon 생성
      - alt amplicon: alt_template_sequence 기반으로 Primer/Amplicon 생성

    핵심:
      - allele-specific primer는 3' 말단을 ref/alt로 강제 치환
      - 치환으로 template에서 exact match가 깨질 수 있으므로 primer3 좌표를 binding_*로 주입
      - reverse primer 3' end 고정은 SEQUENCE_FORCE_RIGHT_START 사용 (중요)
    """

    def __init__(
        self,
        template_sequence: str,
        reference_template_sequence: Optional[str],
        target_start_index: int,
        target_end_index: int,
        target_index: int,
        ref_allele: str,
        alt_allele: str,
        *,
        min_amplicon_length: int = 80,
        max_amplicon_length: int = 100,
        n_primers: int = 50,
        primer3_seq_args: Optional[Dict[str, Any]] = None,
        primer3_global_args: Optional[Dict[str, Any]] = None,
        **kwargs: Any,
    ) -> None:
        self.target_index = int(target_index)
        self.ref_allele = ref_allele.upper()
        self.alt_allele = alt_allele.upper()

        super().__init__(
            template_sequence=template_sequence,  # alt template
            reference_template_sequence=reference_template_sequence,  # ref template
            target_start_index=target_start_index,
            target_end_index=target_end_index,
            min_amplicon_length=min_amplicon_length,
            max_amplicon_length=max_amplicon_length,
            n_primers=n_primers,
            primer3_seq_args=primer3_seq_args,
            primer3_global_args=primer3_global_args,
            **kwargs,
        )

    # ------------------------------------------------------------------
    # primer3 configure
    # ------------------------------------------------------------------
    def _configure_primer_forward_fix(self) -> None:
        super()._configure_primer_common()

        # AS-PCR은 SNP를 primer가 포함해야 하므로 target 회피 규칙 제거
        self.primer3_seq_args.pop("SEQUENCE_TARGET", None)

        # Forward(LEFT) primer 3' end를 target_index에 고정
        self.update_primer3_seq_args({"SEQUENCE_FORCE_LEFT_END": self.target_index})

        self.update_primer3_global_args(
            {
                "PRIMER_PICK_LEFT_PRIMER": 1,
                "PRIMER_PICK_RIGHT_PRIMER": 1,
                "PRIMER_NUM_RETURN": int(self.n_primers),
                "PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
            }
        )

        # 다른 force 제거(혼선 방지)
        self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_START", None)
        self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_END", None)

    def _configure_primer_reverse_fix(self) -> None:
        super()._configure_primer_common()
        self.primer3_seq_args.pop("SEQUENCE_TARGET", None)

        # ✅ RIGHT primer의 3' end를 target_index에 고정하려면 START를 고정해야 함
        self.update_primer3_seq_args({"SEQUENCE_FORCE_RIGHT_END": self.target_index})

        self.update_primer3_global_args(
            {
                "PRIMER_PICK_LEFT_PRIMER": 1,
                "PRIMER_PICK_RIGHT_PRIMER": 1,
                "PRIMER_NUM_RETURN": int(self.n_primers),
                "PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
            }
        )

        # 혼선 방지: 반대 force 제거
        self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_START", None)
        self.primer3_seq_args.pop("SEQUENCE_FORCE_LEFT_END", None)

    # ------------------------------------------------------------------
    # run
    # ------------------------------------------------------------------
    def design(self) -> List[Amplicon]:
        out: List[Amplicon] = []
        # out.extend(self._run_mode_and_build("forward_fix"))
        out.extend(self._run_mode_and_build("reverse_fix"))
        self.amplicon_list = out
        return out

    def _run_mode_and_build(self, mode: str) -> List[Amplicon]:
        self.reset()

        if mode == "forward_fix":
            self._configure_primer_forward_fix()
        elif mode == "reverse_fix":
            self._configure_primer_reverse_fix()
        else:
            raise ValueError(mode)

        self.primer3_result = primer3.bindings.designPrimers(
            seq_args=self.primer3_seq_args,
            global_args=self.primer3_global_args,
        )
        return self._build_amplicons_from_result(self.primer3_result or {}, mode=mode)

    # ------------------------------------------------------------------
    # build amplicons
    # ------------------------------------------------------------------
    def _build_amplicons(self) -> List[Amplicon]:
        if not self.primer3_result:
            return []
        return self._build_amplicons_from_result(self.primer3_result, mode="unknown")

    def _build_amplicons_from_result(self, res: Dict[str, Any], *, mode: str) -> List[Amplicon]:
        n_pairs = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0) or 0)
        if n_pairs == 0:
            return []

        alt_tpl = self.template_sequence            # 변이 반영 template (ALT)
        ref_tpl = self.reference_template_sequence  # ref 반영 template (REF)

        amplicons: List[Amplicon] = []

        for i in range(n_pairs):
            left_seq = res.get(f"PRIMER_LEFT_{i}_SEQUENCE")
            right_seq = res.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
            if not left_seq or not right_seq:
                continue

            # primer3 좌표(0-based within template)
            try:
                left_bind_s, left_bind_e = _pos_to_binding_strict(res.get(f"PRIMER_LEFT_{i}"))
                right_bind_s, right_bind_e = _pos_to_binding_strict(res.get(f"PRIMER_RIGHT_{i}"))
            except Exception:
                continue

            # allele-specific sequences
            fw_ref = _replace_3prime_base(left_seq, self.ref_allele, strand="forward")
            fw_alt = _replace_3prime_base(left_seq, self.alt_allele, strand="forward")

            rv_ref = _replace_3prime_base(right_seq, self.ref_allele, strand="reverse")
            rv_alt = _replace_3prime_base(right_seq, self.alt_allele, strand="reverse")

            # -------------------------
            # forward_fix: forward가 allele-specific, reverse는 공통
            # -------------------------
            if mode == "forward_fix":
                ref_forward_primer = _mk_primer(
                    template_sequence=ref_tpl,
                    reference_template_sequence=ref_tpl,
                    sequence=fw_ref,
                    strand="forward",
                    primer_type="forward",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=left_bind_s,
                    binding_end_index=left_bind_e,
                )
                common_reverse_primer_ref = _mk_primer(
                    template_sequence=ref_tpl,
                    reference_template_sequence=ref_tpl,
                    sequence=right_seq,
                    strand="reverse",
                    primer_type="reverse",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=right_bind_s,
                    binding_end_index=right_bind_e,
                )

                alt_forward_primer = _mk_primer(
                    template_sequence=alt_tpl,
                    reference_template_sequence=ref_tpl,
                    sequence=fw_alt,
                    strand="forward",
                    primer_type="forward",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=left_bind_s,
                    binding_end_index=left_bind_e,
                )
                common_reverse_primer_alt = _mk_primer(
                    template_sequence=alt_tpl,
                    reference_template_sequence=ref_tpl,
                    sequence=right_seq,
                    strand="reverse",
                    primer_type="reverse",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=right_bind_s,
                    binding_end_index=right_bind_e,
                )

                amplicons.append(
                    Amplicon(
                        template_sequence=ref_tpl,
                        reference_template_sequence=ref_tpl,
                        target_start_index=self.target_start_index,
                        target_end_index=self.target_end_index,
                        forward_primer=ref_forward_primer,
                        reverse_primer=common_reverse_primer_ref,
                        assay=f"as_pcr::{mode}",
                        allele="ref",
                    )
                )
                amplicons.append(
                    Amplicon(
                        template_sequence=alt_tpl,
                        reference_template_sequence=ref_tpl,
                        target_start_index=self.target_start_index,
                        target_end_index=self.target_end_index,
                        forward_primer=alt_forward_primer,
                        reverse_primer=common_reverse_primer_alt,
                        assay=f"as_pcr::{mode}",
                        allele="alt",
                    )
                )

            # -------------------------
            # reverse_fix: reverse가 allele-specific, forward는 공통
            # -------------------------
            elif mode == "reverse_fix":
                common_forward_primer_ref = _mk_primer(
                    template_sequence=ref_tpl,
                    reference_template_sequence=ref_tpl,
                    sequence=left_seq,
                    strand="forward",
                    primer_type="forward",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=left_bind_s,
                    binding_end_index=left_bind_e,
                )
                ref_reverse_primer = _mk_primer(
                    template_sequence=ref_tpl,
                    reference_template_sequence=ref_tpl,
                    sequence=rv_ref,
                    strand="reverse",
                    primer_type="reverse",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=right_bind_s,
                    binding_end_index=right_bind_e,
                )

                common_forward_primer_alt = _mk_primer(
                    template_sequence=alt_tpl,
                    reference_template_sequence=ref_tpl,
                    sequence=left_seq,
                    strand="forward",
                    primer_type="forward",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=left_bind_s,
                    binding_end_index=left_bind_e,
                )
                alt_reverse_primer = _mk_primer(
                    template_sequence=alt_tpl,
                    reference_template_sequence=ref_tpl,
                    sequence=rv_alt,
                    strand="reverse",
                    primer_type="reverse",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=right_bind_s,
                    binding_end_index=right_bind_e,
                )

                amplicons.append(
                    Amplicon(
                        template_sequence=ref_tpl,
                        reference_template_sequence=ref_tpl,
                        target_start_index=self.target_start_index,
                        target_end_index=self.target_end_index,
                        forward_primer=common_forward_primer_ref,
                        reverse_primer=ref_reverse_primer,
                        assay=f"as_pcr::{mode}",
                        allele="ref",
                    )
                )
                amplicons.append(
                    Amplicon(
                        template_sequence=alt_tpl,
                        reference_template_sequence=ref_tpl,
                        target_start_index=self.target_start_index,
                        target_end_index=self.target_end_index,
                        forward_primer=common_forward_primer_alt,
                        reverse_primer=alt_reverse_primer,
                        assay=f"as_pcr::{mode}",
                        allele="alt",
                    )
                )

        return amplicons
