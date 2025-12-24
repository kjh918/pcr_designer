
import pysam

from pcr.components.variant import Variant
from pcr.seq.variant import apply_variant

from dataclasses import dataclass
from typing import Optional, Tuple, Dict
import pysam
from Bio.Seq import Seq

# --- 1. 유틸리티 함수: 3' 말단 염기 교체 ---
def _replace_3prime_base(primer_seq: str, allele: str, strand: str) -> str:
    """
    프라이머의 3' 말단을 해당 Allele에 맞게 교체합니다.
    - Forward: Allele 그대로 사용
    - Reverse: Allele의 상보(Complement) 염기 사용 (5'->3' 방향 유지)
    """
    primer_list = list(primer_seq)
    if strand == "forward":
        primer_list[-1] = allele.upper()
    else:
        # Reverse 프라이머는 Sense 가닥의 상보 가닥이므로 상보 염기 필요
        primer_list[-1] = str(Seq(allele).complement()).upper()
    return "".join(primer_list)

# --- 2. 유틸리티 함수: 특이도 강화를 위한 Mismatch 삽입 ---
def _apply_as_pcr_logic(
    primer_seq: str, 
    target_base: str, 
    strand: str, 
    mismatch_pos: int = 3, 
    intensity: str = "strong"
) -> str:
    """
    AS-PCR 프라이머의 특이도를 높이기 위해 3' 말단에서 n번째 위치에 Mismatch를 삽입합니다.
    """
    # 1. 일단 3' 말단을 Allele에 맞게 변경
    primer = _replace_3prime_base(primer_seq, target_base, strand)
    primer_list = list(primer)
    
    # 2. Mismatch 위치 계산 (3' 끝이 1번 위치)
    idx = -mismatch_pos
    original_base = primer_list[idx]
    
    # 3. 강한 Mismatch(Strong) 유도 (Purine <-> Purine, Pyrimidine <-> Pyrimidine 교체 등)
    mismatch_map = {
        "A": "G" if intensity == "strong" else "C",
        "G": "A" if intensity == "strong" else "T",
        "C": "T" if intensity == "strong" else "A",
        "T": "C" if intensity == "strong" else "G"
    }
    
    # 해당 위치의 염기를 강제로 mismatch 염기로 교환
    primer_list[idx] = mismatch_map.get(original_base.upper(), "N")
    return "".join(primer_list)

# --- 3. 유틸리티 함수: 프라이머 서열을 템플릿에 패치 ---
def _patch_template_for_pair(
    template: str, 
    left_seq_5to3: str, left_start: int, left_end: int,
    right_seq_5to3: str, right_start: int, right_end: int
) -> str:
    """
    계산된 프라이머 서열을 템플릿의 바인딩 위치에 덮어씌워(Patch) 
    실제 PCR 반응 시의 '가상 템플릿'을 만듭니다.
    """
    tpl_list = list(template)
    
    # Forward Patch (Same as primer sequence)
    # left_start/end가 0-based 인덱스라고 가정
    tpl_list[left_start : left_end + 1] = list(left_seq_5to3)
    
    # Reverse Patch (Primer의 Reverse Complement를 템플릿에 덮어씌움)
    rc_right = str(Seq(right_seq_5to3).reverse_complement())
    tpl_list[right_start : right_end + 1] = list(rc_right)
    
    return "".join(tpl_list)

# --- 4. 메인 처리 로직 (프라이머 데이터 생성) ---

def process_as_pcr_data(self, left_seq, right_seq, left_bind_s, left_bind_e, right_bind_s, right_bind_e, ref_tpl, alt_tpl):
    # 1. 기본 Allele-specific 프라이머 생성
    fw_ref = _replace_3prime_base(left_seq, self.ref_allele, strand="forward")
    fw_alt = _replace_3prime_base(left_seq, self.alt_allele, strand="forward")

    rv_ref = _replace_3prime_base(right_seq, self.ref_allele, strand="reverse")
    rv_alt = _replace_3prime_base(right_seq, self.alt_allele, strand="reverse")
    
    # 2. Mismatch가 적용된 프라이머 (특이도 강화형)
    fw_ref_mismatch = _apply_as_pcr_logic(left_seq, target_base=self.ref_allele, strand='forward', mismatch_pos=3, intensity="strong") 
    fw_alt_mismatch = _apply_as_pcr_logic(left_seq, target_base=self.alt_allele, strand='forward', mismatch_pos=3, intensity="strong") 
    rv_ref_mismatch = _apply_as_pcr_logic(right_seq, target_base=self.ref_allele, strand='reverse', mismatch_pos=3, intensity="strong") 
    rv_alt_mismatch = _apply_as_pcr_logic(right_seq, target_base=self.alt_allele, strand='reverse', mismatch_pos=3, intensity="strong")

    # 3. 템플릿 패칭 (시뮬레이션용)
    # WT(Ref) 템플릿에 Ref용 프라이머 쌍이 붙은 상태
    wt_tpl = _patch_template_for_pair(
        ref_tpl,
        left_seq_5to3=fw_ref, left_start=left_bind_s, left_end=left_bind_e,
        right_seq_5to3=right_seq, right_start=right_bind_s, right_end=right_bind_e,
    )
    
    # Alt 템플릿에 Alt용 Forward 프라이머가 붙은 상태
    alt_tpl_patched = _patch_template_for_pair(
        alt_tpl,
        left_seq_5to3=fw_alt, left_start=left_bind_s, left_end=left_bind_e,
        right_seq_5to3=right_seq, right_start=right_bind_s, right_end=right_bind_e,
    )
    
    # Alt 템플릿에 Alt용 Reverse 프라이머가 붙은 상태
    rv_alt_tpl_patched = _patch_template_for_pair(
        alt_tpl,
        left_seq_5to3=fw_ref, left_start=left_bind_s, left_end=left_bind_e,
        right_seq_5to3=rv_alt, right_start=right_bind_s, right_end=right_bind_e,
    )

    return {
        "fw_primers": {"ref": fw_ref, "alt": fw_alt, "alt_mm": fw_alt_mismatch},
        "rv_primers": {"ref": rv_ref, "alt": rv_alt, "alt_mm": rv_alt_mismatch},
        "patched_templates": {"wt": wt_tpl, "alt_fw": alt_tpl_patched, "alt_rv": rv_alt_tpl_patched}
    }
    
@dataclass(frozen=True)
class GenomicRegion:
    chrom: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    name: str

# --- 핵심 로직 함수 ---

def _get_contig_length(fasta: pysam.FastaFile, chrom: str) -> Optional[int]:
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
    AS-PCR 용 template window 계산 (1-based inclusive)
    Forward 고정을 위한 Upstream 공간과 Reverse 고정을 위한 Downstream 공간을 모두 확보합니다.
    """
    flank_for_primer = max_primer_length - 1 

    # 1. Forward-fix 시: 타겟 왼쪽(flank) + 타겟 오른쪽(max_amplicon) 필요
    f_up = flank_for_primer
    f_down = (max_amplicon_length + max_primer_length) if (ensure_pair_room and max_amplicon_length) else flank_for_primer
    f_start, f_end = _clamp_window(target_start - f_up, target_end + f_down, contig_length)

    # 2. Reverse-fix 시: 타겟 오른쪽(flank) + 타겟 왼쪽(max_amplicon) 필요
    r_down = flank_for_primer
    r_up = (max_amplicon_length + max_primer_length) if (ensure_pair_room and max_amplicon_length) else flank_for_primer
    r_start, r_end = _clamp_window(target_start - r_up, target_end + r_down, contig_length)

    # 3. 두 경우를 모두 커버하는 UNION 영역 반환
    return min(f_start, r_start), max(f_end, r_end)

def fetch_template_sequence(
    fasta: pysam.FastaFile,
    region: GenomicRegion,
    *,
    max_amplicon_length: Optional[int],
    max_primer_length: int = 30,
    ensure_pair_room: bool = True,
) -> Tuple[str, int, int, int, int]:
    """
    pysam을 이용해 서열을 가져오고 Primer3용 0-based 인덱스를 계산합니다.
    """
    contig_len = _get_contig_length(fasta, region.chrom)

    template_start, template_end = compute_template_window_as_pcr(
        region.start, region.end,
        max_amplicon_length=max_amplicon_length,
        max_primer_length=max_primer_length,
        contig_length=contig_len,
        ensure_pair_room=ensure_pair_room,
    )

    # pysam.fetch(0-based start inclusive, 0-based end exclusive)
    template_sequence = fasta.fetch(
        region.chrom, template_start - 1, template_end
    ).upper()

    # ✅ [중요] 타겟 인덱스 계산 (0-based)
    # 예: region.start가 10이고 template_start가 1이면, 인덱스는 9입니다.
    target_start_index = region.start - template_start
    target_len = region.end - region.start + 1
    target_end_index = target_start_index + target_len - 1

    return template_sequence, template_start, template_end, target_start_index, target_end_index

def build_ref_alt_templates(
    fasta: pysam.FastaFile,
    region: GenomicRegion,
    ref_allele: str,
    alt_allele: str,
    max_amplicon_length: int,
    max_primer_length: int = 30,
    ensure_pair_room: bool = True,
    debug: bool = False,
) -> Dict[str, object]:
    
    (
        reference_sequence,
        t_start, t_end,
        t_start_idx, t_end_idx
    ) = fetch_template_sequence(
        fasta=fasta, region=region,
        max_amplicon_length=max_amplicon_length,
        max_primer_length=max_primer_length,
        ensure_pair_room=ensure_pair_room
    )

    # Variant 적용
    var_ref = Variant(index=t_start_idx, ref=ref_allele, alt=ref_allele)
    var_alt = Variant(index=t_start_idx, ref=ref_allele, alt=alt_allele)

    ref_template_sequence = apply_variant(reference_sequence, var_ref, allele="alt")
    alt_template_sequence = apply_variant(reference_sequence, var_alt, allele="alt")

    if debug:
        print(f"[DEBUG] Window: {region.chrom}:{t_start}-{t_end}")
        print(f"[DEBUG] Target Index in Template: {t_start_idx}")
        print(f"[DEBUG] Base at Index (Ref/Alt): {ref_template_sequence[t_start_idx]}/{alt_template_sequence[t_start_idx]}")

    return {
        "reference_template_sequence": reference_sequence,
        "ref_template_sequence": ref_template_sequence,
        "alt_template_sequence": alt_template_sequence,
        "target_start_index": t_start_idx,
        "target_end_index": t_end_idx,
        "template_start": t_start,
        "template_end": t_end,
    }