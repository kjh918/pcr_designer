from __future__ import annotations
from typing import Literal

Meth = Literal["M", "U"]

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
