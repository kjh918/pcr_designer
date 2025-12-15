from __future__ import annotations

from typing import Any, Dict, List, Optional, Literal

from designers.base import BasePrimerDesigner
from models import Primer, Probe, Amplicon, Variant


Allele = Literal["ref", "alt"]
Meth = Literal["M", "U"]  # MSP-style bisulfite model


def apply_variant(reference_seq: str, var: Variant, *, allele: Allele) -> str:
    s = reference_seq.upper()
    if s[var.index] != var.ref.upper():
        raise ValueError(
            f"Reference base mismatch at {var.index}: expected {var.ref.upper()}, got {s[var.index]}"
        )
    base = var.ref.upper() if allele == "ref" else var.alt.upper()
    return s[:var.index] + base + s[var.index + 1:]


def bisulfite_convert(seq: str, *, meth: Meth) -> str:
    """
    - M: CpG의 C 유지, non-CpG C는 T
    - U: 모든 C를 T
    """
    s = seq.upper()
    out = list(s)
    for i, ch in enumerate(out):
        if ch != "C":
            continue
        is_cpg = (i + 1 < len(out) and out[i + 1] == "G")
        if meth == "M" and is_cpg:
            out[i] = "C"
        else:
            out[i] = "T"
    return "".join(out)


class BisulfiteVariantProbeDesigner(BasePrimerDesigner):
    """
    Variant assay + bisulfite + (internal oligo) probe.
    - reference_template_sequence: 항상 ref allele 기준(reference_gdna) 고정
    - template_sequence: allele(ref/alt) 적용 후 bisulfite 변환된 템플릿
    - probe는 variant 1개를 반드시 덮고 allele-specific base를 가져야 함
    """

    def __init__(
        self,
        reference_gdna: str,  # ref allele 기준 gDNA (reference_template_sequence로 고정)
        variant: Variant,
        target_start_index: int,
        target_end_index: int,
        *,
        allele: Allele = "ref",
        meth: Meth = "U",
        # probe 조건(primer3 internal oligo)
        pick_probe: bool = True,
        primer3_seq_args: Optional[Dict[str, Any]] = None,
        primer3_global_args: Optional[Dict[str, Any]] = None,
        **kwargs,
    ) -> None:
        self.variant = variant
        self.allele = allele
        self.meth = meth
        self.assay = "bisulfite_variant_probe"

        # (중요) reference는 항상 ref allele 고정
        reference_template_sequence = reference_gdna

        # template은 allele만 바꿔서 생성
        allele_gdna = apply_variant(reference_gdna, variant, allele=allele)

        # bisulfite 변환
        template_sequence = bisulfite_convert(allele_gdna, meth=meth)

        super().__init__(
            template_sequence=template_sequence,
            reference_template_sequence=reference_template_sequence,
            target_start_index=target_start_index,
            target_end_index=target_end_index,
            primer3_seq_args=primer3_seq_args,
            primer3_global_args=primer3_global_args,
            **kwargs,
        )

        if pick_probe:
            self.update_primer3_global_args({"PRIMER_PICK_INTERNAL_OLIGO": 1})

    def _build_amplicons(self) -> List[Amplicon]:
        if not self.primer3_result:
            return []

        res = self.primer3_result
        n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        out: List[Amplicon] = []

        for i in range(n):
            left_seq = res.get(f"PRIMER_LEFT_{i}_SEQUENCE")
            right_seq = res.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
            if not left_seq or not right_seq:
                continue

            left_start, left_len = res[f"PRIMER_LEFT_{i}"]
            right_start, right_len = res[f"PRIMER_RIGHT_{i}"]

            fwd = Primer(
                template_sequence=self.template_sequence,
                reference_template_sequence=self.reference_template_sequence,  # ref allele 고정
                sequence=left_seq,
                strand="forward",
                primer_type="forward",
                target_start_index=self.target_start_index,
                target_end_index=self.target_end_index,
                start_index=left_start,
                length=left_len,
                assay=self.assay,
                allele=self.allele,
                converted=True,
            )
            rev = Primer(
                template_sequence=self.template_sequence,
                reference_template_sequence=self.reference_template_sequence,
                sequence=right_seq,
                strand="reverse",
                primer_type="reverse",
                target_start_index=self.target_start_index,
                target_end_index=self.target_end_index,
                start_index=right_start,
                length=right_len,
                assay=self.assay,
                allele=self.allele,
                converted=True,
            )

            probe_obj = None
            probe_seq = res.get(f"PRIMER_INTERNAL_{i}_SEQUENCE")
            if probe_seq:
                p_start, p_len = res[f"PRIMER_INTERNAL_{i}"]
                probe_obj = Probe(
                    template_sequence=self.template_sequence,
                    reference_template_sequence=self.reference_template_sequence,
                    sequence=probe_seq,
                    strand="forward",  # internal oligo는 template 방향으로 해석(필요 시 확장)
                    primer_type="probe",
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    start_index=p_start,
                    length=p_len,
                    assay=self.assay,
                    allele=self.allele,
                    converted=True,
                    variant=self.variant,
                )

                # 핵심: probe는 반드시 변이를 덮고 allele-specific 이어야 함
                if not probe_obj.validate_variant_specific():
                    continue

            out.append(
                Amplicon(
                    template_sequence=self.template_sequence,
                    reference_template_sequence=self.reference_template_sequence,
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    forward_primer=fwd,
                    reverse_primer=rev,
                    probe=probe_obj,
                    assay=self.assay,
                    allele=self.allele,
                )
            )

        return out


def design_ref_alt(
    reference_gdna: str,
    variant: Variant,
    target_start_index: int,
    target_end_index: int,
    **kwargs,
) -> Dict[str, List[Amplicon]]:
    """
    ref/alt를 각각 template_sequence만 바꿔서 디자인 결과를 반환
    """
    ref_designer = BisulfiteVariantProbeDesigner(
        reference_gdna,
        variant,
        target_start_index,
        target_end_index,
        allele="ref",
        **kwargs,
    )
    alt_designer = BisulfiteVariantProbeDesigner(
        reference_gdna,
        variant,
        target_start_index,
        target_end_index,
        allele="alt",
        **kwargs,
    )
    return {"ref": ref_designer.design(), "alt": alt_designer.design()}
