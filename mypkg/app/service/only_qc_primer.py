
from primer.qc import SimpleAmplicon
from primer.qc import QCThresholds, evaluate_amplicons, make_amplicon_for_qc
from primer.pcr_components import Primer, Amplicon  # 방금 올린 Primer/Amplicon 정의

def make_amplicon_for_qc(
    forward_seq: str,
    reverse_seq: str,
    probe_seq: str | None = None,
    template_seq: str | None = None,
) -> Amplicon:
    """
    QC only용 Amplicon 생성:
    - template_sequence가 주어지면 그걸 사용
    - 없으면 forward + 'N'*10 + reverse 를 이어붙인 가짜 템플릿 생성
    """
    forward_seq = forward_seq.strip().upper()
    reverse_seq = reverse_seq.strip().upper()
    probe_seq = probe_seq.strip().upper() if probe_seq else None

    if template_seq and template_seq.strip():
        template = template_seq.strip().upper()
    else:
        # QC-only 용 가짜 템플릿 (둘 다 포함되도록)
        template = forward_seq + ("N" * 10) + reverse_seq

    target_start_index = 0
    target_end_index = len(template) - 1

    # Primer 객체 생성 (Primer.__init__ 안에서 Tm, GC, hairpin, homodimer 다 계산됨)
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
            strand="forward",
            primer_type="probe",
            target_start_index=target_start_index,
            target_end_index=target_end_index,
        )

    amp = Amplicon(
        template_sequence=template,
        reference_template_sequence=template,
        target_start_index=target_start_index,
        target_end_index=target_end_index,
        forward_primer=f_primer,
        reverse_primer=r_primer,
        probe=probe_primer,
    )

    return amp

def evaluate_primer_set_with_components(
    genomic_id: str,
    forward_seq: str,
    reverse_seq: str,
    probe_seq: str | None,
    template_seq: str | None,
    qc_thresholds: QCThresholds,
) -> dict:
    amp = make_amplicon_for_qc(
        forward_seq=forward_seq,
        reverse_seq=reverse_seq,
        probe_seq=probe_seq,
        template_seq=template_seq,
    )

    total_rows, filtered_rows = evaluate_amplicons(
        genomic_id=genomic_id,
        amplicons=[amp],
        qc_thresholds=qc_thresholds,
    )
    return total_rows[0]  # 한 세트만 넘겼으니까 첫 번째 row 사용
