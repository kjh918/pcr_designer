#!/usr/bin/env python3
import sys
import os
import json
import argparse
from typing import Dict, Any, List, Optional

try:
    import pysam
except ImportError:
    pysam = None

# 상위 경로 추가 (pcr 패키지 인식용)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import load_pipeline_config
from pcr.factory import PCRFactory
from pcr.designers.base.schema import BaseDesignInput 

def reverse_complement(seq: str) -> str:
    """염기서열의 역상보 서열을 반환합니다."""
    return seq.translate(str.maketrans('ATGCatgcNn', 'TACGtacgNn'))[::-1]

def fetch_and_validate_ref(chrom: str, start: int, end: int, fasta_path: str, expected_ref: str, strand: str = "+"):
    """FASTA에서 서열을 추출하여 검증합니다."""
    if not pysam:
        raise ImportError("pysam is required (pip install pysam)")
    
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    with pysam.FastaFile(fasta_path) as fasta:
        actual_ref_plus = fasta.fetch(chrom, start - 1, end).upper()
        actual_ref_for_cmp = reverse_complement(actual_ref_plus) if strand == "-" else actual_ref_plus
        
        if actual_ref_for_cmp != expected_ref.upper():
            raise ValueError(
                f"Reference mismatch! FASTA has '{actual_ref_for_cmp}' on '{strand}' strand, "
                f"but input was '{expected_ref}'"
            )
            
    return actual_ref_plus

def _generate_aspcr_templates(
    left_seq: str, right_seq: str, ref: str, alt: str, strand: str, 
    fixed_prime: str = "reverse", 
    mismatch_pos: Optional[int] = 3,
    mismatch_intensity: str = "strong"
) -> Dict[str, str]:
    """방향, 위치, 강도(intensity)에 맞춰 프라이머 서열 길이는 유지한 채 해당 염기만 '치환'합니다."""
    ref_strand = reverse_complement(ref) if strand == "-" else ref
    alt_strand = reverse_complement(alt) if strand == "-" else alt
    
    wt_mm_left, alt_mm_left = left_seq, left_seq
    wt_mm_right, alt_mm_right = right_seq, right_seq
    
    def _apply_as_pcr_logic(target_base: str, intensity: str) -> str:
        mismatch_map = {
            "A": "G" if intensity == "strong" else "C",
            "G": "A" if intensity == "strong" else "T",
            "C": "T" if intensity == "strong" else "A",
            "T": "C" if intensity == "strong" else "G"
        }
        return mismatch_map.get(target_base.upper(), "N")

    # 미스매치 위치가 유효할 때만 치환 (길이는 절대 변하지 않음)
    if mismatch_pos and mismatch_pos > 1:
        if fixed_prime == "forward":
            mm_index = -(mismatch_pos - 1)
            if len(left_seq) >= abs(mm_index):
                original_base = left_seq[mm_index]
                new_base = _apply_as_pcr_logic(original_base, mismatch_intensity)
                wt_mm_left = left_seq[:mm_index] + new_base + left_seq[mm_index+1:]
                alt_mm_left = wt_mm_left
        else:
            mm_index = mismatch_pos - 2
            if len(right_seq) > mm_index:
                original_base = right_seq[mm_index]
                new_base = _apply_as_pcr_logic(original_base, mismatch_intensity)
                wt_mm_right = right_seq[:mm_index] + new_base + right_seq[mm_index+1:]
                alt_mm_right = wt_mm_right

    return {
        "wt": left_seq + ref_strand + right_seq,
        "alt": left_seq + alt_strand + right_seq,
        "wt_mm": wt_mm_left + ref_strand + wt_mm_right,
        "alt_mm": alt_mm_left + alt_strand + alt_mm_right
    }

def design_aspcr_primers(
    chrom: str, start: int, end: int, ref: str, alt: str, strand: str,
    genome: str = "hg38",
    fasta_path: str = None,
    padding: int = 100, 
    top_k: int = 5,
    base_yaml: str = "pcr/config/base_pcr.yaml", 
    system_yaml: str = "pcr/config/system.yaml",
    fixed_prime: str = "reverse", 
    mismatch_pos: int = 3,
    mismatch_intensity: str = "strong",
    qc_overrides: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    
    assay_type = "ASPCR"
    task_name = f"{chrom}_{start}_{ref}>{alt}({strand})"

    # 1. Config 로드
    user_overrides = {}
    if qc_overrides:
        user_overrides["qc_criteria"] = qc_overrides

    config = load_pipeline_config(base_yaml, system_yaml, assay_type, user_overrides=user_overrides)

    # 2. FASTA 경로 결정
    if not fasta_path:
        try:
            ref_cfg = config.get_reference(genome)
            fasta_path = ref_cfg.fasta_path
        except ValueError as e:
            return {"status": "error", "reason": f"Genome config error: {str(e)}"}
        except Exception as e:
            return {"status": "error", "reason": f"Failed to resolve fasta path: {str(e)}"}
    
    if not fasta_path:
        return {"status": "error", "reason": f"FASTA path is required. Check system.yaml for '{genome}'."}

    # 3. 서열 추출 및 템플릿 생성
    try:
        ref_plus = fetch_and_validate_ref(chrom, start, end, fasta_path, ref, strand)
        alt_plus = reverse_complement(alt) if strand == "-" else alt
        
        template_start_0based = max(0, start - 1 - padding)
        template_end_0based = end + padding
        
        with pysam.FastaFile(fasta_path) as fasta:
            left_seq = fasta.fetch(chrom, template_start_0based, start - 1).upper()
            right_seq = fasta.fetch(chrom, end, template_end_0based).upper()
            
        # 🔥 사용자 설정 (위치, 강도) 반영하여 템플릿 생성
        templates = _generate_aspcr_templates(
            left_seq, right_seq, ref, alt, strand, 
            fixed_prime=fixed_prime,
            mismatch_pos=mismatch_pos,
            mismatch_intensity=mismatch_intensity
        )
        rel_start = len(left_seq)
        rel_end = rel_start + len(alt)
        
    except Exception as e:
        return {"status": "error", "reason": f"Sequence extraction failed: {str(e)}"}

    # 4. Factory 실행
    factory = PCRFactory(config)
    output = factory.run(
        assay_type=assay_type,
        name=task_name,
        target_start=rel_start,
        target_end=rel_end,
        reference_name=genome,
        template_sequence=templates["alt"],
        reference_sequence=templates["wt"],
        top_k=top_k,
        run_qc=True,
        overrides={"PRIMER_NUM_RETURN": 5},
        template_genomic_start=template_start_0based, # 절대좌표 시작점 명시적 전달
        target_strand=strand,
        templates=templates,
        fixed_prime=fixed_prime,
        mismatch_pos=mismatch_pos if mismatch_pos else None 
    )

    if output.status != "success":
        fail_reason = output.error_msg or "Failed to design AS-PCR primers."
        logs = getattr(output, "log_messages", [])
        if logs:
            log_str = "; ".join(logs) if isinstance(logs, list) else str(logs)
            if len(log_str) > 500: log_str = log_str[:500] + "..."
            return {"status": "fail", "reason": fail_reason, "log_messages": log_str}
        return {"status": "fail", "reason": fail_reason}

    # 5. 결과 포맷팅
    def format_region(oligo):
        if not oligo or getattr(oligo, "region", None) is None:
            return {"genomic_pos": "Unknown", "strand": "+", "index_start": None}
        reg = oligo.region
        return {
            "genomic_pos": f"{reg.chrom}:{reg.start}-{reg.end}",
            "strand": reg.strand,
            "index_start": oligo.start_index,
            "index_end": oligo.start_index + len(oligo.sequence)
        }

    sets_dict = {}
    for amp in output.amplicons: 
        set_id = getattr(amp, "set_id", "UnknownSet")
        if set_id not in sets_dict:
            sets_dict[set_id] = []
        sets_dict[set_id].append(amp)

    results = []
    rank = 1

    # 🔥 수정된 포맷팅: 덮어쓰기 방지 및 세트 레벨 래핑 적용
    for set_id, amplicons_in_set in sets_dict.items():
        set_qc_pass = all(getattr(amp, "is_qc_pass", False) for amp in amplicons_in_set) and len(amplicons_in_set) == 4
        rep_amp = amplicons_in_set[0]
        alignment_data = getattr(rep_amp, "alignment_visual", ["Alignment data not available."])
        
        set_data = {
            "rank": rank,
            "set_id": set_id,
            "set_qc_pass": set_qc_pass, 
            "fixed_prime": fixed_prime,
            "common_metrics": {
                "pair_penalty": round(rep_amp.pair_penalty, 3),
                "product_size": rep_amp.product_size if hasattr(rep_amp, 'product_size') else None
            },
            "alignment_text_block": "\n".join(alignment_data) if alignment_data else "",
            "alleles": {} # 여기에 wt, alt, wt_mm, alt_mm 담김
        }
        
        for amp in amplicons_in_set:
            # 💡 Pydantic 제약 우회: id에서 allele_type 파싱
            if "wt_mm" in amp.id:
                allele_type = "wt_mm"
            elif "alt_mm" in amp.id:
                allele_type = "alt_mm"
            elif "wt" in amp.id:
                allele_type = "wt"
            else:
                allele_type = "alt"
            
            fwd = amp.forward
            rev = amp.reverse
            
            fwd_tm = getattr(fwd, "tm", 0.0)
            fwd_gc = getattr(fwd, "gc_percent", getattr(fwd, "gc", 0.0))
            fwd_penalty = getattr(fwd, "penalty", 0.0)
            
            rev_tm = getattr(rev, "tm", 0.0)
            rev_gc = getattr(rev, "gc_percent", getattr(rev, "gc", 0.0))
            rev_penalty = getattr(rev, "penalty", 0.0)
            
            pair_penalty = getattr(amp, "pair_penalty", fwd_penalty + rev_penalty + abs(fwd_tm - rev_tm))
            amp_tm = getattr(amp, "tm", 0.0)
            amp_gc = getattr(amp, "gc_percent", getattr(amp, "gc", 0.0))

            amp_seq = amp.template_sequence[fwd.start_index : rev.end_index] if amp.template_sequence else ""

            set_data["alleles"][allele_type] = {
                "id": amp.id,
                "pair_penalty": round(pair_penalty, 2),
                "qc_info": {
                    "is_pass": getattr(amp, "is_qc_pass", False), 
                    "fail_reason": getattr(amp, "qc_log", getattr(amp, "qc_fail_reason", "None")),
                    "blast_details": getattr(amp, "blast_stats", None)
                },
                "amplicon_info": {
                    "sequence": amp_seq,
                    "length": len(amp_seq),
                    "tm": round(amp_tm, 2) if amp_tm else None, 
                    "gc": round(amp_gc, 2) if amp_gc else None,
                    "genomic_pos": f"{chrom}:{fwd.region.start}-{rev.region.end}" if getattr(fwd, "region", None) else "Unknown"
                },
                "oligos": {
                    "forward": {
                        "sequence": fwd.sequence,
                        "tm": round(fwd_tm, 2),
                        "gc": round(fwd_gc, 2),
                        "penalty": round(fwd_penalty, 2),
                        "hairpin_tm": round(getattr(fwd, "hairpin_tm", 0.0), 2),     
                        "homodimer_tm": round(getattr(fwd, "homodimer_tm", 0.0), 2), 
                        "is_allele_specific": getattr(fwd, "is_allele_specific", False),
                        "terminal_base": getattr(fwd, "terminal_base", None),
                        "mismatch_base": getattr(fwd, "mismatch_base", None),
                        **format_region(fwd)
                    },
                    "reverse": {
                        "sequence": rev.sequence,
                        "tm": round(rev_tm, 2),
                        "gc": round(rev_gc, 2),
                        "penalty": round(rev_penalty, 2),
                        "hairpin_tm": round(getattr(rev, "hairpin_tm", 0.0), 2),
                        "homodimer_tm": round(getattr(rev, "homodimer_tm", 0.0), 2),
                        "is_allele_specific": getattr(rev, "is_allele_specific", False),
                        "terminal_base": getattr(rev, "terminal_base", None),
                        "mismatch_base": getattr(rev, "mismatch_base", None),
                        **format_region(rev)
                    }
                }
            }
        
        results.append(set_data)
        rank += 1

    return {
        "status": "success",
        "metadata": {
            "assay": "AS-PCR", "chrom": chrom, "valid_sets_found": len(results),
            "genome_build": genome, "fasta_used": fasta_path
        },
        "target_info": {
            "input_strand": strand, "input_ref": ref, "input_alt": alt,
            "index_start": rel_start, "index_end": rel_end,
            "genomic_pos": f"{chrom}:{start}-{end}",
        },
        "results": results
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Design AS-PCR primers")
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--alt", required=True)
    parser.add_argument("--strand", choices=["+", "-"], default="+")
    parser.add_argument("--genome", default="hg38")
    parser.add_argument("--fasta", help="Optional override")
    parser.add_argument("-p", "--padding", type=int, default=150)
    parser.add_argument("-k", "--top_k", type=int, default=5)
    parser.add_argument("--base_config", default="pcr/config/base_pcr.yaml")
    parser.add_argument("--system_config", default="pcr/config/system.yaml")
    
    # 🔥 옵션들 추가 완료
    parser.add_argument("--fixed_prime", choices=["forward", "reverse"], default="forward", help="Anchor side for the SNP")
    parser.add_argument("-m", "--mismatch_pos", type=int, default=3, help="Mismatch position from 3' end (e.g. 2, 3). Use 0 for none.")
    parser.add_argument("--intensity", choices=["strong", "weak"], default="strong", help="Mismatch intensity (strong/weak)")
    
    args = parser.parse_args()
    
    res = design_aspcr_primers(
        chrom=args.chrom, start=args.start, end=args.end, ref=args.ref, alt=args.alt, strand=args.strand,
        genome=args.genome, fasta_path=args.fasta, padding=args.padding, top_k=args.top_k,
        base_yaml=args.base_config, system_yaml=args.system_config,
        fixed_prime=args.fixed_prime,
        mismatch_pos=args.mismatch_pos,       
        mismatch_intensity=args.intensity     
    )
    
    print(json.dumps(res, indent=2))