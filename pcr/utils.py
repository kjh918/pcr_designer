
from typing import Tuple
from Bio.Seq import reverse_complement

def get_start_end_index(template_sequence: str, sequence: str) -> Tuple[int, int]:
    """
    template_sequence에서 sequence(또는 reverse complement)의 시작/끝 인덱스를 반환 (0-based, end inclusive)

    NOTE: 반복 서열이 있으면 index()는 첫 매치로 잡힘.
          primer3가 준 start/length가 있으면 그걸 우선 쓰는 설계를 권장.
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
