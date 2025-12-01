# primer/core.py

from typing import Optional

import pysam

from app.schemas import RegionInput, PrimerPair, Primer
from primer.qpcr_designer import qPCRdesigner
from primer.qpcr_designer import qPCRdesigner

from primer.reference import get_fasta_handle   # 혹은 _get_fasta_handle

# 만약 config.settings 에 reference path가 들어있다면:
#   settings.references["hg19"].fasta 이런 식으로 접근한다고 가정
try:
    from config.settings import settings
except ImportError:
    settings = None  # 필요하면 나중에 직접 경로 넘겨도 됨


def _primer_obj_to_schema(primer_obj) -> Optional[Primer]:
    """
    pcr_components.Primer / Amplicon 내 primer 객체를
    app.schemas.Primer 로 변환.
    """
    if primer_obj is None:
        return None

    if hasattr(primer_obj, "sequence"):
        seq = primer_obj.sequence
    elif hasattr(primer_obj, "seq"):
        seq = primer_obj.seq
    else:
        raise ValueError("Primer object must have 'sequence' or 'seq' attribute")

    tm = getattr(primer_obj, "tm", None)
    gc = getattr(primer_obj, "gc", None)
    start = getattr(primer_obj, "start", None)
    end = getattr(primer_obj, "end", None)
    strand = getattr(primer_obj, "strand", None)

    return Primer(
        seq=seq,
        tm=tm,
        gc=gc,
        start=start,
        end=end,
        strand=strand,
    )



# ============================================================
#   qPCR 전용 wrapper : qPCRdesigner 사용
# ============================================================

def design_qpcr_for_region(
    region: RegionInput,
    reference_name: str,
    min_amplicon_length: int = 80,
    max_amplicon_length: int = 120,
    n_probes: int = 100,
    n_primers: int = 100,
    bisulfite: bool = False,
    n_best: int = 10,           # ★ 화면에 보여줄 후보 개수
):
    fasta = get_fasta_handle(reference_name)

    qp = qPCRdesigner(
        f_reference_fasta=fasta,
        chrom=region.chrom,
        start=region.start,
        end=region.end,
        min_amplicon_length=min_amplicon_length,
        max_amplicon_length=max_amplicon_length,
        n_probes=n_probes,
        n_primers=n_primers,
        bisulfite=bisulfite,
    )

    qp.design_primer()
    
    if not qp.amplicon_list:
        raise ValueError("qPCRdesigner가 어떤 amplicon도 생성하지 못했습니다.")

    # 1) 최상위 1개 → 기존처럼 PrimerPair 로
    best_amp = qp.amplicon_list[0]

    forward = _primer_obj_to_schema(getattr(best_amp, "forward_primer", None))
    reverse = _primer_obj_to_schema(getattr(best_amp, "reverse_primer", None))
    probe   = _primer_obj_to_schema(getattr(best_amp, "probe", None))

    product_size = getattr(best_amp, "amplicon_length", None)
    if product_size is None:
        product_size = getattr(best_amp, "product_size", None)

    best_pair = PrimerPair(
        forward=forward,
        reverse=reverse,
        probe=probe,
        product_size=product_size,
    )

    # 2) 상위 n_best 개를 테이블용 row 리스트로
    primer_rows = []
    for rank, amp in enumerate(qp.amplicon_list[:n_best], start=1):
        f = getattr(amp, "forward_primer", None)
        r = getattr(amp, "reverse_primer", None)
        p = getattr(amp, "probe", None)

        row = {
            "rank": rank,
            "product_size": getattr(amp, "amplicon_length", None) or getattr(amp, "product_size", None),

            "forward_seq": getattr(f, "sequence", None),
            "forward_tm": getattr(f, "tm", None),
            "forward_gc": getattr(f, "gc", None),

            "reverse_seq": getattr(r, "sequence", None),
            "reverse_tm": getattr(r, "tm", None),
            "reverse_gc": getattr(r, "gc", None),

            "probe_seq": getattr(p, "sequence", None),
            "probe_tm": getattr(p, "tm", None),
            "probe_gc": getattr(p, "gc", None),

            # 선택 기준들 있으면 getattr 로 추가 (없으면 None)
            "probe_cpg_count": getattr(amp, "probe_cpg_count", None),
            "forward_cpg_count": getattr(amp, "forward_cpg_count", None),
            "reverse_cpg_count": getattr(amp, "reverse_cpg_count", None),
        }
        primer_rows.append(row)

    # ★ main 에서는 (최상위 1쌍, N개 row 리스트)를 같이 받도록
    return best_pair, primer_rows