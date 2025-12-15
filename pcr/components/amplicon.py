from __future__ import annotations

from typing import Any, Dict, Optional, Literal

from pcr.components.primer import Primer
from pcr.utils import get_start_end_index

Allele = Literal["ref", "alt"]

class Amplicon:
    def __init__(
        self,
        template_sequence: str,
        target_start_index: int,
        target_end_index: int,
        reference_template_sequence: Optional[str] = None,
        chrom: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        forward_primer: Optional[Primer] = None,
        reverse_primer: Optional[Primer] = None,
        probe: Optional[Primer] = None,
        assay: str = "generic",
        allele: Allele = "ref",
    ) -> None:
        self.template_sequence = template_sequence
        self.reference_template_sequence = reference_template_sequence or template_sequence

        self.target_start_index = target_start_index
        self.target_end_index = target_end_index

        self.chrom = chrom
        self.start = start
        self.end = end

        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.probe = probe

        self.assay = assay
        self.allele = allele

        if self.forward_primer is not None:
            self.forward_start_index, self.forward_end_index = get_start_end_index(
                self.template_sequence, self.forward_primer.sequence
            )
        if self.reverse_primer is not None:
            self.reverse_start_index, self.reverse_end_index = get_start_end_index(
                self.template_sequence, self.reverse_primer.sequence
            )
        if self.probe is not None:
            self.probe_start_index, self.probe_end_index = get_start_end_index(
                self.template_sequence, self.probe.sequence
            )

        self.amplicon_sequence = self._calc_amplicon_sequence()

    def _calc_amplicon_sequence(self) -> Optional[str]:
        if self.forward_primer is None or self.reverse_primer is None:
            return None
        return self.template_sequence[self.forward_start_index : self.reverse_end_index + 1]

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "reference_template_sequence": self.reference_template_sequence,
            "template_sequence": self.template_sequence,
            "target_start_index": self.target_start_index,
            "target_end_index": self.target_end_index,
            "assay": self.assay,
            "allele": self.allele,
        }

        if self.amplicon_sequence is not None:
            d["amplicon_sequence"] = self.amplicon_sequence
            d["amplicon_length"] = len(self.amplicon_sequence)
        else:
            d["amplicon_sequence"] = None
            d["amplicon_length"] = None

        if self.forward_primer is not None:
            d.update(self.forward_primer.to_dict())
        if self.reverse_primer is not None:
            d.update(self.reverse_primer.to_dict())
        if self.probe is not None:
            d.update(self.probe.to_dict())
        return d
