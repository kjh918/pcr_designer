# primer/qc/factories.py
from __future__ import annotations

from typing import Optional

from pcr.components import Primer, Amplicon


def make_amplicon_for_qc(
    forward_seq: str,
    reverse_seq: str,
    probe_seq: str | None = None,
    template_seq: str | None = None,
) -> Amplicon:
    forward_seq = forward_seq.replace(' ','').strip().upper()
    reverse_seq = reverse_seq.replace(' ','').strip().upper()
    probe_seq = probe_seq.replace(' ','').strip().upper() if probe_seq else None

    if template_seq and template_seq.strip():
        template = template_seq.strip().upper().replace(' ','')
    else:
        template = forward_seq + ("N" * 10) + reverse_seq

    if probe_seq is not None:
        template = forward_seq + ("N" * 10) + probe_seq + ("N" * 10) + reverse_seq

    target_start_index = 0
    target_end_index = len(template) - 1

    f_primer = Primer(
        template_sequence=template,
        reference_template_sequence=template,
        sequence=forward_seq,
        strand="forward",
        primer_type="forward",
        target_start_index=target_start_index,
        target_end_index=target_end_index,
    )

    r_primer = Primer(
        template_sequence=template,
        reference_template_sequence=template,
        sequence=reverse_seq,
        strand="reverse",
        primer_type="reverse",
        target_start_index=target_start_index,
        target_end_index=target_end_index,
    )

    probe_primer: Primer | None = None
    if probe_seq:
        probe_primer = Primer(
            template_sequence=template,
            reference_template_sequence=template,
            sequence=probe_seq,
            strand="reverse",
            primer_type="probe",
            target_start_index=target_start_index,
            target_end_index=target_end_index,
        )

    return Amplicon(
        template_sequence=template,
        reference_template_sequence=template,
        target_start_index=target_start_index,
        target_end_index=target_end_index,
        forward_primer=f_primer,
        reverse_primer=r_primer,
        probe=probe_primer,
    )
