# primer/reference.py
from typing import Tuple
from pyfaidx import Fasta
import pysam
from config.settings import settings  # 참고: 여기서 fasta 경로 불러온다고 가정

# ref_name(hg19/hg38) 별 Fasta 핸들 캐시
_FASTA_CACHE: dict[str, Fasta] = {}


def get_fasta_handle(ref_name: str) -> Fasta:
    """
    config.settings.references[ref_name].fasta 를 이용해서
    pyfaidx.Fasta 객체를 가져온다 (캐시 사용).
    """
    if ref_name not in _FASTA_CACHE:
        ref_cfg = settings.references.get(ref_name)
        if ref_cfg is None:
            raise ValueError(f"Unknown reference name: {ref_name}")
        _FASTA_CACHE[ref_name] = pysam.FastaFile(str(ref_cfg.fasta))
    return _FASTA_CACHE[ref_name]


def fetch_sequence_from_reference(
    ref_name: str,
    chrom: str,
    start: int,
    end: int,
    max_amplicon_length: int,
) -> Tuple[str, int]:
    """
    qPCRdesigner 스타일로 reference 에서 서열을 가져온다.

    - 입력:
        ref_name: 'hg19', 'hg38' 등
        chrom: 'chr1', '1' 등 (FASTA와 일치해야 함)
        start, end: 타겟 구간의 게놈 좌표 (1-based, start <= end)
        max_amplicon_length: 좌우로 펼칠 길이 (qPCRdesigner의 max_amplicon_length)

    - 동작:
        실제 fetch 구간 = [end - max_amplicon_length, start + max_amplicon_length + 1)
        (pyfaidx에서는 end가 exclusive)

    - 출력:
        (template_sequence, template_start)
        template_sequence: 펼쳐진 템플릿 서열 (대문자)
        template_start: 템플릿 첫 번째 base의 게놈 좌표 (1-based)
    """

    if start is None or end is None:
        raise ValueError("start / end 좌표가 필요합니다.")
    if start > end:
        raise ValueError(f"start({start}) > end({end}) 인 잘못된 영역입니다.")

    fasta = get_fasta_handle(ref_name)

    # 크롬존 존재 여부 확인
    try:
        chrom_seq = fasta[chrom]
    except KeyError:
        raise ValueError(f"Reference '{ref_name}' 에 chrom '{chrom}' 이(가) 없습니다.")

    chrom_len = len(chrom_seq)

    # qPCRdesigner 패턴 그대로:
    raw_fetch_start = end - max_amplicon_length
    raw_fetch_end = start + max_amplicon_length + 1  # pyfaidx end exclusive

    # 영역이 통째로 reference 바깥이면 에러
    if raw_fetch_end < 1 or raw_fetch_start > chrom_len:
        raise ValueError(
            f"요청한 영역이 reference 범위를 완전히 벗어났습니다: "
            f"fetch_start={raw_fetch_start}, fetch_end={raw_fetch_end}, chr_len={chrom_len}"
        )

    # reference 안으로 클램핑
    fetch_start = max(raw_fetch_start, 1)
    fetch_end = min(raw_fetch_end, chrom_len)

    # pyfaidx: [start-1 : end] (end는 exclusive)
    template_seq = str(chrom_seq[fetch_start - 1:fetch_end]).upper()
    template_start = fetch_start  # 템플릿 첫 base의 게놈 좌표
    template_end = fetch_end

    # 안전장치: 템플릿 길이가 너무 짧으면 경고/에러
    if len(template_seq) < (end - start + 1):
        raise ValueError(
            f"템플릿 길이({len(template_seq)})가 타겟 길이({end - start + 1})보다 짧습니다. "
            f"입력 좌표 또는 max_amplicon_length를 확인하세요."
        )

    return template_seq, template_start, template_end
