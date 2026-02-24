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

def design_qpcr_primers(
    chrom: str, start: int, end: int, ref: str, alt: str, strand: str,
    genome: str = "hg38",
    fasta_path: str = None,
    padding: int = 100, 
    top_k: int = 5,
    base_yaml: str = "pcr/config/base_pcr.yaml", 
    system_yaml: str = "pcr/config/system.yaml",
    qc_overrides: Optional[Dict[str, Any]] = None # (Optional) 추가됨
) -> Dict[str, Any]:
    
    assay_type = "qpcr"
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

    # 3. 서열 추출
    try:
        ref_plus = fetch_and_validate_ref(chrom, start, end, fasta_path, ref, strand)
        alt_plus = reverse_complement(alt) if strand == "-" else alt
        
        template_start_0based = max(0, start - 1 - padding)
        
        with pysam.FastaFile(fasta_path) as fasta:
            left_seq = fasta.fetch(chrom, template_start_0based, start - 1).upper()
            right_seq = fasta.fetch(chrom, end, end + padding).upper()
            
        alt_template = left_seq + alt_plus.upper() + right_seq
        ref_template = left_seq + ref_plus.upper() + right_seq
        
        rel_start = len(left_seq)
        rel_end = rel_start + len(alt_plus)
        
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
        template_sequence=alt_template,
        reference_sequence=ref_template, 
        top_k=top_k,
        run_qc=True,
        overrides={"PRIMER_NUM_RETURN": 10},
        template_genomic_start=template_start_0based,
        target_strand=strand
    )
    print(output)
    # 🚨 [개선] 실패 시 로그 메시지까지 포함해서 반환
    if output.status != "success":
        fail_reason = output.error_msg or "Failed to design primers/probes."
        
        # 로그가 있다면 합쳐서 디버깅을 돕습니다.
        logs = getattr(output, "log_messages", [])
        if logs:
            # 리스트면 join, 문자열이면 그대로
            log_str = "; ".join(logs) if isinstance(logs, list) else str(logs)
            # 너무 길면 자름
            if len(log_str) > 500: log_str = log_str[:500] + "..."
            
            # 최종 반환값에 log_messages 필드 추가 (main.py에서 출력 가능)
            return {
                "status": "fail", 
                "reason": fail_reason,
                "log_messages": log_str
            }
            
        return {
            "status": "fail", 
            "reason": fail_reason
        }

    # 5. 결과 포맷팅
    def get_position_info(seq: str, is_reverse: bool = False, is_probe: bool = False):
        if not seq: return None
        search_seq = reverse_complement(seq) if is_reverse else seq
        rel_idx = alt_template.find(search_seq)
        strand_sign = "+" if not is_reverse else "-"
        if is_probe and rel_idx == -1:
            search_seq = reverse_complement(seq)
            rel_idx = alt_template.find(search_seq)
            strand_sign = "-" if rel_idx != -1 else "+"
        if rel_idx == -1:
            return {"index_start": None, "index_end": None, "genomic_pos": "Unknown", "strand": strand_sign}
        length = len(seq)
        genomic_start = template_start_0based + rel_idx + 1 
        genomic_end = genomic_start + length - 1
        return {
            "index_start": rel_idx,
            "index_end": rel_idx + length,
            "strand": strand_sign,
            "genomic_pos": f"{chrom}:{genomic_start}-{genomic_end}"
        }

    results = []
    for rank, amp in enumerate(output.amplicons, start=1):
        alignment_data = getattr(amp, "alignment_visual", ["Alignment data not available."])
        fwd_pos = get_position_info(amp.forward.sequence, is_reverse=False)
        rev_pos = get_position_info(amp.reverse.sequence, is_reverse=True)
        prb_pos = get_position_info(amp.probe.sequence if amp.probe else "", is_probe=True)

        results.append({
            "rank": rank,
            "id": amp.id,
            "metrics": {"pair_penalty": round(amp.pair_penalty, 3)},
            "qc_info": {"is_pass": amp.is_qc_pass, "fail_reason": getattr(amp, "qc_fail_reason", "None")},
            "amplicon_info": {
                "sequence": amp.sequence, "length": amp.product_size,
                "tm": round(amp.tm, 2), "gc": amp.gc, "genomic_pos": amp.genomic_pos
            },
            "oligos": {
                "forward": {
                    "sequence": amp.forward.sequence, "length": len(amp.forward.sequence),
                    "tm": round(amp.forward.tm, 2), "gc": round((amp.forward.sequence.count('G') + amp.forward.sequence.count('C')) / len(amp.forward.sequence) * 100, 2),
                    **fwd_pos
                },
                "reverse": {
                    "sequence": amp.reverse.sequence, "length": len(amp.reverse.sequence),
                    "tm": round(amp.reverse.tm, 2), "gc": round((amp.reverse.sequence.count('G') + amp.reverse.sequence.count('C')) / len(amp.reverse.sequence) * 100, 2),
                    **rev_pos
                },
                "probe": {
                    "sequence": amp.probe.sequence, "length": len(amp.probe.sequence),
                    "tm": round(amp.probe.tm, 2), "gc": round((amp.probe.sequence.count('G') + amp.probe.sequence.count('C')) / len(amp.probe.sequence) * 100, 2),
                    **prb_pos
                } if amp.probe else None
            },
            "alignment_text_block": "\n".join(alignment_data)
        })

    return {
        "status": "success",
        "metadata": {
            "assay": "TaqMan-qPCR", "chrom": chrom, "candidates_found": len(results),
            "genome_build": genome, "fasta_used": fasta_path
        },
        "target_info": {
            "input_strand": strand, "input_ref": ref, "input_alt": alt,
            "mapped_plus_strand_ref": ref_plus, "mapped_plus_strand_alt": alt_plus,
            "index_start": rel_start, "index_end": rel_end,
            "genomic_pos": f"{chrom}:{start}-{end}", "template_sequence": alt_template 
        },
        "results": results
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Design qPCR primers/probes")
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--alt", required=True)
    parser.add_argument("--strand", choices=["+", "-"], default="+")
    parser.add_argument("--genome", default="hg38")
    parser.add_argument("--fasta", help="Optional override")
    parser.add_argument("-p", "--padding", type=int, default=100)
    parser.add_argument("-k", "--top_k", type=int, default=5)
    parser.add_argument("--base_config", default="pcr/config/base_pcr.yaml")
    parser.add_argument("--system_config", default="pcr/config/system.yaml")
    
    args = parser.parse_args()
    
    res = design_qpcr_primers(
        chrom=args.chrom, start=args.start, end=args.end, ref=args.ref, alt=args.alt, strand=args.strand,
        genome=args.genome, fasta_path=args.fasta, padding=args.padding, top_k=args.top_k,
        base_yaml=args.base_config, system_yaml=args.system_config
    )
    
    print(json.dumps(res, indent=2))