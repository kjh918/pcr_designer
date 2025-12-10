
from primer.qc import SimpleAmplicon
from primer.qc import QCThresholds, evaluate_amplicons, make_amplicon_for_qc
from primer.pcr_components import Primer, Amplicon  # 방금 올린 Primer/Amplicon 정의


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
    return total_rows[0]  #
