from __future__ import annotations
from typing import Literal
from pcr.components.variant import Variant

Allele = Literal["ref", "alt"]

def apply_variant(reference_seq: str, variant: Variant, *, allele: Allele) -> str:
    """
    reference_seq: ref allele 기준 서열 (reference_template_sequence로 고정할 대상)
    allele: "ref" or "alt"
    """
    s = reference_seq.upper()
    if s[variant.index] != variant.ref.upper():
        raise ValueError(
            f"Reference base mismatch at {variant.index}: expected {variant.ref.upper()}, got {s[variant.index]}"
        )
    base = variant.ref.upper() if allele == "ref" else variant.alt.upper()
    print(s[:variant.index] + base + s[variant.index + 1:])
    print(base)
    
    return s[:variant.index] + base + s[variant.index + 1:]
