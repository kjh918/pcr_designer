from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Tuple
import pysam

@dataclass(frozen=True)
class GenomicRegion:
    chrom: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive

def compute_template_window(
    target_start: int,
    target_end: int,
    *,
    max_amplicon_length: Optional[int],
) -> Tuple[int, int]:
    """target(1-based inclusive) -> template window(1-based inclusive)"""
    target_len = target_end - target_start + 1

    if max_amplicon_length is not None and max_amplicon_length > target_len:
        extra = max_amplicon_length - target_len
        upstream_extra = extra // 2
        downstream_extra = extra - upstream_extra
    else:
        upstream_extra = 0
        downstream_extra = 0

    template_start = max(1, target_start - upstream_extra)
    template_end = target_end + downstream_extra
    return template_start, template_end

def fetch_template_sequence(
    fasta: pysam.FastaFile,
    region: GenomicRegion,
    *,
    max_amplicon_length: Optional[int],
) -> Tuple[str, int, int, int, int]:
    """
    Returns:
      template_sequence (upper)
      template_start (1-based inclusive)
      template_end   (1-based inclusive)
      target_start_index (0-based within template)
      target_end_index   (0-based within template)
    """
    template_start, template_end = compute_template_window(
        region.start, region.end, max_amplicon_length=max_amplicon_length
    )

    template_sequence = fasta.fetch(
        region.chrom,
        template_start - 1,  # 0-based
        template_end,        # end param is treated as end-exclusive in pysam
    ).upper()

    target_len = region.end - region.start + 1
    target_start_index = region.start - template_start
    target_end_index = target_start_index + target_len - 1

    return template_sequence, template_start, template_end, target_start_index, target_end_index
