from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, Dict

import pysam

from pcr.components.variant import Variant
from pcr.seq.variant import apply_variant


@dataclass(frozen=True)
class GenomicRegion:
    chrom: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    name: str


def _get_contig_length(fasta: pysam.FastaFile, chrom: str) -> Optional[int]:
    """contig 길이를 알 수 있으면 반환. 모르면 None."""
    try:
        return int(fasta.get_reference_length(chrom))
    except Exception:
        return None


def _clamp_window(start: int, end: int, contig_length: Optional[int]) -> Tuple[int, int]:
    start = max(1, start)
    if contig_length is not None:
        end = min(contig_length, end)
    if end < start:
        end = start
    return start, end


def compute_template_window_as_pcr(
    target_start: int,
    target_end: int,
    *,
    max_amplicon_length: Optional[int],
    max_primer_length: int = 30,
    contig_length: Optional[int] = None,
    ensure_pair_room: bool = True,
) -> Tuple[int, int]:
    """
    AS-PCR 용 template window (1-based inclusive)

    - forward 고정 설계 가능하도록 upstream 확보
    - reverse 고정 설계 가능하도록 downstream 확보
    - pair(product size)까지 고려해 한쪽에 충분한 room 확보(옵션)
    - forward-fix window + reverse-fix window를 합집합으로 결합(UNION)
    """
    if max_primer_length <= 0:
        raise ValueError("max_primer_length must be > 0")

    flank_for_primer = max_primer_length - 1  # primer 3' end 고정 시 최소 flank

    # forward-fix에 필요한 window
    f_up = flank_for_primer
    if ensure_pair_room and max_amplicon_length is not None:
        f_down = max_amplicon_length + max_primer_length
    else:
        f_down = flank_for_primer

    f_start, f_end = _clamp_window(target_start - f_up, target_end + f_down, contig_length)

    # reverse-fix에 필요한 window
    r_down = flank_for_primer
    if ensure_pair_room and max_amplicon_length is not None:
        r_up = max_amplicon_length + max_primer_length
    else:
        r_up = flank_for_primer

    r_start, r_end = _clamp_window(target_start - r_up, target_end + r_down, contig_length)

    # UNION
    template_start, template_end = _clamp_window(min(f_start, r_start), max(f_end, r_end), contig_length)
    return template_start, template_end


def fetch_template_sequence(
    fasta: pysam.FastaFile,
    region: GenomicRegion,
    *,
    bisulfite: bool = False,
    max_amplicon_length: Optional[int],
    max_primer_length: int = 30,
    ensure_pair_room: bool = True,
) -> Tuple[str, int, int, int, int]:
    """
    Returns:
      template_sequence (upper)
      template_start (1-based inclusive)
      template_end   (1-based inclusive)
      target_start_index (0-based within template)
      target_end_index   (0-based within template)

    ✅ AS-PCR forward/reverse 고정 모두 커버하도록 compute_template_window_as_pcr() 사용
    """
    contig_len = _get_contig_length(fasta, region.chrom)

    template_start, template_end = compute_template_window_as_pcr(
        region.start,
        region.end,
        max_amplicon_length=max_amplicon_length,
        max_primer_length=max_primer_length,
        contig_length=contig_len,
        ensure_pair_room=ensure_pair_room,
    )

    # pysam.fetch: start=0-based inclusive, end=0-based exclusive
    template_sequence = fasta.fetch(
        region.chrom,
        template_start - 1,
        template_end,
    ).upper()
    # target index within template (0-based)
    target_len = region.end - region.start + 1
    target_start_index = region.start - template_start
    target_end_index = target_start_index + target_len - 1
    print(template_sequence[target_start_index])


    # bisulfite는 여기서는 미적용(필요하면 후처리 훅 추가)
    _ = bisulfite

    return template_sequence, template_start , template_end, target_start_index, target_end_index


def build_ref_alt_templates(
    *,
    fasta: pysam.FastaFile,
    region: GenomicRegion,
    ref_allele: str,
    alt_allele: str,
    max_amplicon_length: int,
    max_primer_length: int = 30,
    ensure_pair_room: bool = True,
    debug: bool = False,
) -> Dict[str, object]:
    """
    Returns dict:
      - reference_template_sequence
      - ref_template_sequence
      - alt_template_sequence
      - target_start_index / target_end_index (0-based in template)
      - template_start / template_end (1-based genomic window)
    """
    (
        reference_template_sequence,
        template_start,
        template_end,
        target_start_index,
        target_end_index,
    ) = fetch_template_sequence(
        fasta=fasta,
        region=region,
        max_amplicon_length=max_amplicon_length,
        max_primer_length=max_primer_length,
        ensure_pair_room=ensure_pair_room,
    )

    # Variant 정의 (template 기준 index)
    variant = Variant(
        index=target_start_index,
        ref=ref_allele,
        alt=alt_allele,
    )
    print(variant)
    ref_template_sequence = apply_variant(reference_template_sequence, variant, allele="ref")
    alt_template_sequence = apply_variant(reference_template_sequence, variant, allele="alt")
    
    # ✅ 서열 변화 확인 (타겟 위치 염기 확인)
    if debug:
        i = target_start_index
        ref_base_expected = ref_allele.upper()
        alt_base_expected = alt_allele.upper()

        ref_base_actual = ref_template_sequence[i].upper()
        alt_base_actual = alt_template_sequence[i].upper()
        base_ref0 = reference_template_sequence[i].upper()

        print("[DEBUG] template window:", region.chrom, template_start, template_end, "len", len(reference_template_sequence))
        print("[DEBUG] target index (0-based in template):", i)
        print("[DEBUG] reference_template base:", base_ref0)
        print("[DEBUG] ref_template base:", ref_base_actual, "expected", ref_base_expected)
        print("[DEBUG] alt_template base:", alt_base_actual, "expected", alt_base_expected)

        if ref_base_actual != ref_base_expected:
            print("[WARN] ref_template base mismatch at target_index")
        if alt_base_actual != alt_base_expected:
            print("[WARN] alt_template base mismatch at target_index")

        # “진짜로 바뀌었는지”도 체크 (SNP면 보통 ref!=alt)
        if ref_base_actual == alt_base_actual:
            print("[WARN] ref/alt templates have same base at target_index (check variant/ref-alt input)")

    return dict(
        reference_template_sequence=reference_template_sequence,
        ref_template_sequence=ref_template_sequence,
        alt_template_sequence=alt_template_sequence,
        target_start_index=target_start_index,
        target_end_index=target_end_index,
        template_start=template_start,
        template_end=template_end,
    )
