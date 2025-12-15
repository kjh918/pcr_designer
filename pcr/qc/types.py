# primer/qc/types.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Any, Optional, TypedDict


@dataclass
class SimpleAmplicon:
    forward_sequence: str
    reverse_sequence: str
    probe_sequence: str | None = None

    def to_dict(self) -> Dict[str, Any]:
        # 구현은 factories/thermo에서 주입할 수도 있지만,
        # 기존 구조를 최대한 유지하려면 thermo.evaluate_amplicons에서
        # compute_* 호출해서 채우는 쪽이 더 깔끔함.
        return {
            "forward_sequence": self.forward_sequence,
            "reverse_sequence": self.reverse_sequence,
            "probe_sequence": self.probe_sequence,
        }


class BlastHit(TypedDict):
    qseqid: str
    sseqid: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qseq: str
    sseq: str
