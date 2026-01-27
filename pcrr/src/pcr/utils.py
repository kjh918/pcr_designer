from typing import Tuple, Optional
from Bio.Seq import reverse_complement

def get_start_end_index(template_sequence: str, sequence: str) -> Tuple[int, int]:
    """
    반복 서열 대응:
    - 정방향은 첫 매치
    - 역방향(RC)은 마지막 매치 사용
    """
    tpl = template_sequence.upper()
    seq = sequence.upper()

    # 1) forward 그대로
    idx = tpl.find(seq)
    if idx != -1:
        return idx, idx + len(seq) - 1

    # 2) reverse primer → RC를 마지막 매치로
    rc = reverse_complement(seq).upper()
    idx = tpl.rfind(rc)   # ⭐ 핵심 변경
    if idx != -1:
        return idx, idx + len(seq) - 1

    raise ValueError(f"{sequence} or its reverse complement not found in template.")
