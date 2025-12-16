# components/as_assay.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Any, Optional

from pcr.components.primer import Primer

@dataclass(frozen=True)
class AsPcrAssay:
    wt_forward: Primer
    alt_forward: Primer
    reverse: Primer
    amplicon_start_index: int
    amplicon_end_index: int
    amplicon_sequence: str

    def to_dict(self) -> Dict[str, Any]:
        d = {}
        d.update(self.wt_forward.to_dict())
        d.update(self.alt_forward.to_dict())
        d.update(self.reverse.to_dict())
        d["amplicon_start_index"] = self.amplicon_start_index
        d["amplicon_end_index"] = self.amplicon_end_index
        d["amplicon_length"] = len(self.amplicon_sequence)
        d["amplicon_sequence"] = self.amplicon_sequence
        return d
