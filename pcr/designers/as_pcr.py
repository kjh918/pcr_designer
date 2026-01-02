from __future__ import annotations

from typing import Any, Dict, List, Tuple
import primer3
from Bio.Seq import Seq

from pcr.designers.base import BasePrimerDesigner
from pcr.components import Amplicon
from pcr.components.primer import Primer

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def _replace_3prime_base(primer_5to3: str, base: str, strand: str = "forward") -> str:
    s = list(primer_5to3.upper())
    if not s:
        return ""
    if strand == "forward":
        s[-1] = base.upper()
    else:
        s[-1] = str(Seq(base).complement()).upper()
    return "".join(s)


def _apply_as_pcr_logic(
    primer_5to3: str,
    target_base: str,
    strand: str,
    mismatch_pos: int = 3,
    intensity: str = "strong",
) -> str:
    mismatch_map = {
        "strong": {"A": "G", "T": "C", "G": "A", "C": "T"},
        "medium": {"A": "T", "T": "A", "G": "C", "C": "G"},
        "weak": {"A": "C", "T": "G", "G": "T", "C": "A"},
    }

    primer_seq = _replace_3prime_base(primer_5to3, target_base, strand)
    s = list(primer_seq)

    idx = -int(mismatch_pos)
    if abs(idx) <= len(s):
        original = s[idx]
        new_base = mismatch_map.get(intensity.lower(), mismatch_map["strong"]).get(original, "A")
        if new_base == original:
            for cand in ("A", "C", "G", "T"):
                if cand != original:
                    new_base = cand
                    break
        s[idx] = new_base

    return "".join(s)


def _patch_template_for_pair(
    template_5to3: str,
    *,
    left_seq_5to3: str,
    left_start: int,
    left_end: int,
    right_seq_5to3: str,
    right_start: int,
    right_end: int,
) -> str:
    """Overwrites the template with actual primer sequences (Right is RC'd)."""
    t = list(template_5to3.upper())
    left_seq = left_seq_5to3.upper()
    right_on_template = str(Seq(right_seq_5to3.upper()).reverse_complement())

    t[left_start : left_end + 1] = list(left_seq)
    t[right_start : right_end + 1] = list(right_on_template)
    return "".join(t)


def _mk_primer(
    *,
    template_sequence: str,
    reference_template_sequence: str,
    sequence: str,
    strand: str,
    primer_type: str,
    target_start_index: int,
    target_end_index: int,
    binding_start_index: int,
    binding_end_index: int,
) -> Primer:
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
    if not pos or len(pos) != 2:
        raise ValueError(f"Invalid position: {pos}")
    start, length = int(pos[0]), int(pos[1])
    return start, start + length - 1


# ----------------------------------------------------------------------
# MODIFIED: AS-PCR "3' 고정 + 길이 가변" 후보를 직접 생성 (Right primer용)
# ----------------------------------------------------------------------

def _enumerate_right_candidates_3prime_fixed(
    template_5to3: str,
    *,
    target_index: int,
    min_len: int,
    max_len: int,
) -> List[Tuple[str, Tuple[int, int]]]:
    """
    Right primer(5'->3') candidates where the primer's 3' end binds to template[target_index].
    Implementation:
      - Take template window [target_index : target_index + L] (forward strand)
      - Right primer sequence is reverse-complement of that window
      - Binding on template is (start=target_index, end=target_index+L-1)
    """
    t = template_5to3.upper()
    out: List[Tuple[str, Tuple[int, int]]] = []
    for L in range(min_len, max_len + 1):
        start = target_index
        end = target_index + L - 1
        if start < 0 or end >= len(t):
            continue
        window = t[start : end + 1]
        seq = str(Seq(window).reverse_complement())  # primer 5'->3'
        out.append((seq, (start, end)))
    return out


def _gc_percent(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100.0


# ----------------------------------------------------------------------
# AsPcrDesigner
# ----------------------------------------------------------------------

class AsPcrDesigner(BasePrimerDesigner):
    # primer3-py built-in max primer length = 36
    P3_MAX_PRIMER_LEN = 36  # MODIFIED

    def __init__(
        self,
        template_sequence: str,
        reference_template_sequence: str,
        target_start_index: int,
        target_end_index: int,
        target_index: int,
        ref_allele: str,
        alt_allele: str,
        **kwargs,
    ) -> None:
        self.target_index = int(target_index)
        self.ref_allele = ref_allele.upper()
        self.alt_allele = alt_allele.upper()
        super().__init__(
            template_sequence=template_sequence,
            reference_template_sequence=reference_template_sequence,
            target_start_index=target_start_index,
            target_end_index=target_end_index,
            **kwargs,
        )

    def reset(self) -> None:
        self.primer3_seq_args = {"SEQUENCE_TEMPLATE": self.template_sequence}
        self.primer3_global_args = {
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_MAX_SIZE": 25,
            "PRIMER_OPT_TM": 55.0,
            "PRIMER_MIN_TM": 35.0,
            "PRIMER_MAX_TM": 95.0,
            "PRIMER_OPT_GC": 45.0,
            "PRIMER_MIN_GC": 20.0,
            "PRIMER_MAX_GC": 85.0,
            "PRIMER_MAX_POLY_X": 5,
            "PRIMER_SALT_MONOVALENT": 50.0,
            "PRIMER_DNA_CONC": 50.0,
            "PRIMER_EXPLAIN_FLAG": 1,
        }

    def _configure_primer_forward_fix(self) -> None:
        super()._configure_primer_common()
        self.primer3_seq_args.pop("SEQUENCE_TARGET", None)

        # LEFT primer 3' end fixed (OK)
        self.update_primer3_seq_args({"SEQUENCE_FORCE_LEFT_END": int(self.target_index)})

        self.update_primer3_global_args(
            {
                "PRIMER_PICK_LEFT_PRIMER": 1,
                "PRIMER_PICK_RIGHT_PRIMER": 1,
                "PRIMER_NUM_RETURN": int(self.n_primers),
                "PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
            }
        )

        self.primer3_seq_args.pop("SEQUENCE_FORCE_LEFT_START", None)
        self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_START", None)
        self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_END", None)

    # ------------------------------------------------------------------
    # MODIFIED: reverse_fix는 primer3 FORCE_RIGHT_*를 믿지 않고,
    #           "3' 고정 + 길이 가변" right primer 후보를 직접 만들어
    #           primer3에는 그 right를 고정(SEQUENCE_PRIMER_REVCOMP)으로 넣고 left만 설계하게 함
    # ------------------------------------------------------------------
    def _design_reverse_fix_by_enumeration(self) -> List[Amplicon]:
        # common 세팅
        super()._configure_primer_common()
        self.primer3_seq_args.pop("SEQUENCE_TARGET", None)

        # product size range/num return 유지
        self.update_primer3_global_args(
            {
                "PRIMER_PICK_LEFT_PRIMER": 1,
                "PRIMER_PICK_RIGHT_PRIMER": 0,  # MODIFIED: right는 우리가 고정
                "PRIMER_NUM_RETURN": int(self.n_primers),
                "PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
            }
        )

        # FORCE_* 제거 (혼선 방지)
        self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_START", None)  # MODIFIED
        self.primer3_seq_args.pop("SEQUENCE_FORCE_RIGHT_END", None)
        self.primer3_seq_args.pop("SEQUENCE_FORCE_LEFT_START", None)
        self.primer3_seq_args.pop("SEQUENCE_FORCE_LEFT_END", None)

        # candidate right primers (3' fixed)
        min_len = int(self.primer3_global_args.get("PRIMER_MIN_SIZE", 18))
        max_len = min(int(self.primer3_global_args.get("PRIMER_MAX_SIZE", 25)), self.P3_MAX_PRIMER_LEN)

        cands = _enumerate_right_candidates_3prime_fixed(
            self.template_sequence,
            target_index=self.target_index,
            min_len=min_len,
            max_len=max_len,
        )

        # 디버그: 후보 길이/GC 요약
        print(f"[reverse_fix][DEBUG] Enumerated right candidates count={len(cands)} (len {min_len}-{max_len})")  # MODIFIED
        if cands:
            gc_list = [_gc_percent(s) for s, _ in cands]
            print(f"[reverse_fix][DEBUG] right GC% range: {min(gc_list):.1f} - {max(gc_list):.1f}")  # MODIFIED

        amplicons: List[Amplicon] = []
        ref_tpl_raw = self.reference_template_sequence
        alt_tpl_raw = self.template_sequence

        # 각 right 후보를 고정해서 left를 설계
        for idx, (right_seq, (r_s, r_e)) in enumerate(cands):
            seq_args = dict(self.primer3_seq_args)
            seq_args["SEQUENCE_PRIMER_REVCOMP"] = right_seq  # MODIFIED: right 고정

            res = primer3.bindings.designPrimers(seq_args, self.primer3_global_args) or {}

            # explain 출력(필요 시)
            if idx == 0:
                print("[reverse_fix][Primer3 Explain(sample)]",
                      res.get("PRIMER_LEFT_EXPLAIN"),
                      res.get("PRIMER_RIGHT_EXPLAIN"),
                      res.get("PRIMER_PAIR_EXPLAIN"))

            n_pairs = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))
            if n_pairs == 0:
                continue

            # 여기서 res의 RIGHT는 고정프라이머라서 PRIMER_RIGHT_0_*가 없을 수도 있음.
            # LEFT만 primer3 결과로 받고, right는 우리가 만든 seq/pos를 사용.
            for i in range(n_pairs):
                l_seq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
                l_s, l_e = _pos_to_binding_strict(res[f"PRIMER_LEFT_{i}"])

                # AS-PCR variant 생성 (reverse가 allele-specific)
                rv_ref = _replace_3prime_base(right_seq, self.ref_allele, "reverse")
                rv_alt = _replace_3prime_base(right_seq, self.alt_allele, "reverse")
                rv_ref_mm = _apply_as_pcr_logic(right_seq, self.ref_allele, "reverse", mismatch_pos=3)
                rv_alt_mm = _apply_as_pcr_logic(right_seq, self.alt_allele, "reverse", mismatch_pos=3)

                configs = [
                    ("wt", l_seq, rv_ref, ref_tpl_raw, "ref"),
                    ("alt", l_seq, rv_alt, alt_tpl_raw, "alt"),
                    ("wt_mismatch", l_seq, rv_ref_mm, ref_tpl_raw, "ref"),
                    ("alt_mismatch", l_seq, rv_alt_mm, alt_tpl_raw, "alt"),
                ]

                for label, f_seq, r_seq_final, base_tpl, allele_type in configs:
                    patched_tpl = _patch_template_for_pair(
                        base_tpl,
                        left_seq_5to3=f_seq,
                        left_start=l_s,
                        left_end=l_e,
                        right_seq_5to3=r_seq_final,
                        right_start=r_s,
                        right_end=r_e,
                    )

                    f_primer = _mk_primer(
                        template_sequence=patched_tpl,
                        reference_template_sequence=ref_tpl_raw,
                        sequence=f_seq,
                        strand="forward",
                        primer_type="forward",
                        target_start_index=self.target_start_index,
                        target_end_index=self.target_end_index,
                        binding_start_index=l_s,
                        binding_end_index=l_e,
                    )
                    r_primer = _mk_primer(
                        template_sequence=patched_tpl,
                        reference_template_sequence=ref_tpl_raw,
                        sequence=r_seq_final,
                        strand="reverse",
                        primer_type="reverse",
                        target_start_index=self.target_start_index,
                        target_end_index=self.target_end_index,
                        binding_start_index=r_s,
                        binding_end_index=r_e,
                    )

                    amplicons.append(
                        Amplicon(
                            template_sequence=patched_tpl,
                            reference_template_sequence=ref_tpl_raw,
                            target_start_index=self.target_start_index,
                            target_end_index=self.target_end_index,
                            forward_primer=f_primer,
                            reverse_primer=r_primer,
                            assay=f"as_pcr::reverse_fix::cand{idx}::set{i}::{label}",  # MODIFIED
                            allele=allele_type,
                        )
                    )

        return amplicons

    # ------------------------------------------------------------------

    def design(self) -> List[Amplicon]:
        out: List[Amplicon] = []
        out.extend(self._run_mode_and_build("forward_fix"))
        # MODIFIED: reverse_fix는 enumeration 방식으로 수행
        self.reset()
        out.extend(self._design_reverse_fix_by_enumeration())  # MODIFIED
        self.amplicon_list = out
        return out

    def _run_mode_and_build(self, mode: str) -> List[Amplicon]:
        self.reset()

        if mode == "forward_fix":
            self._configure_primer_forward_fix()
        else:
            # MODIFIED: reverse_fix는 여기로 안 들어오게 처리
            return []

        res = primer3.bindings.designPrimers(self.primer3_seq_args, self.primer3_global_args) or {}
        return self._build_amplicons_from_result(res or {}, mode=mode)

    def _build_amplicons_from_result(self, res: Dict[str, Any], mode: str) -> List[Amplicon]:
        n_pairs = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))
        if n_pairs == 0:
            return []

        amplicons: List[Amplicon] = []
        ref_tpl_raw = self.reference_template_sequence
        alt_tpl_raw = self.template_sequence

        for i in range(n_pairs):
            l_seq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
            r_seq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
            l_s, l_e = _pos_to_binding_strict(res[f"PRIMER_LEFT_{i}"])
            r_s, r_e = _pos_to_binding_strict(res[f"PRIMER_RIGHT_{i}"])

            fw_ref = _replace_3prime_base(l_seq, self.ref_allele, "forward")
            fw_alt = _replace_3prime_base(l_seq, self.alt_allele, "forward")
            fw_ref_mm = _apply_as_pcr_logic(l_seq, self.ref_allele, "forward", mismatch_pos=3)
            fw_alt_mm = _apply_as_pcr_logic(l_seq, self.alt_allele, "forward", mismatch_pos=3)

            if mode == "forward_fix":
                configs = [
                    ("wt", fw_ref, r_seq, ref_tpl_raw, "ref"),
                    ("alt", fw_alt, r_seq, alt_tpl_raw, "alt"),
                    ("wt_mismatch", fw_ref_mm, r_seq, ref_tpl_raw, "ref"),
                    ("alt_mismatch", fw_alt_mm, r_seq, alt_tpl_raw, "alt"),
                ]
            else:
                configs = []

            for label, f_seq, r_seq_final, base_tpl, allele_type in configs:
                patched_tpl = _patch_template_for_pair(
                    base_tpl,
                    left_seq_5to3=f_seq,
                    left_start=l_s,
                    left_end=l_e,
                    right_seq_5to3=r_seq_final,
                    right_start=r_s,
                    right_end=r_e,
                )

                f_primer = _mk_primer(
                    template_sequence=patched_tpl,
                    reference_template_sequence=ref_tpl_raw,
                    sequence=f_seq,
                    strand="forward",
                    primer_type="forward",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=l_s,
                    binding_end_index=l_e,
                )
                r_primer = _mk_primer(
                    template_sequence=patched_tpl,
                    reference_template_sequence=ref_tpl_raw,
                    sequence=r_seq_final,
                    strand="reverse",
                    primer_type="reverse",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    binding_start_index=r_s,
                    binding_end_index=r_e,
                )

                amplicons.append(
                    Amplicon(
                        template_sequence=patched_tpl,
                        reference_template_sequence=ref_tpl_raw,
                        target_start_index=self.target_start_index,
                        target_end_index=self.target_end_index,
                        forward_primer=f_primer,
                        reverse_primer=r_primer,
                        assay=f"as_pcr::{mode}::set{i}::{label}",
                        allele=allele_type,
                    )
                )

        return amplicons
