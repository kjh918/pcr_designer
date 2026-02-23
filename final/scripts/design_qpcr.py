#!/usr/bin/env python3
import sys
import os
import json
import argparse
from typing import Dict, Any, List

try:
    import pysam
except ImportError:
    pysam = None

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import load_pipeline_config
from pcr.factory import PCRFactory

def reverse_complement(seq: str) -> str:
    """염기서열의 역상보 서열을 반환합니다."""
    return seq.translate(str.maketrans('ATGCatgcNn', 'TACGtacgNn'))[::-1]

def fetch_and_validate_ref(chrom: str, start: int, end: int, fasta_path: str, expected_ref: str, strand: str = "+"):
    """FASTA에서 서열을 추출하여 검증합니다."""
    if not pysam:
        raise ImportError("pysam이 필요합니다. (pip install pysam)")
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
    fasta_path: str, padding: int, top_k: int,
    base_yaml: str, system_yaml: str
) -> Dict[str, Any]:
    
    assay_type = "qpcr"
    task_name = f"{chrom}_{start}_{ref}>{alt}({strand})"

    config = load_pipeline_config(base_yaml, system_yaml, assay_type)

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
        return {"status": "error", "reason": str(e)}

    factory = PCRFactory(config)
    output = factory.run(
        assay_type=assay_type,
        name=task_name,
        target_start=rel_start,
        target_end=rel_end,
        reference_name=chrom, # 좌표 계산을 위해 크로모좀 이름 전달
        template_sequence=alt_template,
        reference_sequence=ref_template, 
        top_k=top_k,
        run_qc=True,
        overrides={"PRIMER_NUM_RETURN": 10},
        template_genomic_start=template_start_0based, # 절대 좌표 오프셋 전달
        target_strand=strand
    )

    if output.status != "success":
        return {
            "status": "fail", 
            "reason": output.error_msg or "Failed to design primers/probes."
        }

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
            "metrics": {
                "pair_penalty": round(amp.pair_penalty, 3)
            },
            "qc_info": {
                "is_pass": amp.is_qc_pass,
                "fail_reason": getattr(amp, "qc_fail_reason", "None") 
            },
            "amplicon_info": {
                "sequence": amp.sequence,         
                "length": amp.product_size,       
                "tm": round(amp.tm, 2),           
                "gc": amp.gc,                     
                "genomic_pos": amp.genomic_pos    
            },
            "oligos": {
                "forward": {
                    "sequence": amp.forward.sequence,
                    "length": len(amp.forward.sequence),
                    "tm": round(amp.forward.tm, 2),
                    "gc": round((amp.forward.sequence.count('G') + amp.forward.sequence.count('C')) / len(amp.forward.sequence) * 100, 1),
                    **fwd_pos
                },
                "reverse": {
                    "sequence": amp.reverse.sequence,
                    "length": len(amp.reverse.sequence),
                    "tm": round(amp.reverse.tm, 2),
                    "gc": round((amp.reverse.sequence.count('G') + amp.reverse.sequence.count('C')) / len(amp.reverse.sequence) * 100, 1),
                    **rev_pos
                },
                "probe": {
                    "sequence": amp.probe.sequence,
                    "length": len(amp.probe.sequence),
                    "tm": round(amp.probe.tm, 2),
                    "gc": round((amp.probe.sequence.count('G') + amp.probe.sequence.count('C')) / len(amp.probe.sequence) * 100, 1),
                    **prb_pos
                } if amp.probe else None
            },
            "alignment_text_block": "\n".join(alignment_data)
        })

    return {
        "status": "success",
        "metadata": {
            "assay": "TaqMan-qPCR (SNP)",
            "chrom": chrom,
            "candidates_found": len(results)
        },
        "target_info": {
            "input_strand": strand,
            "input_ref": ref,
            "input_alt": alt,
            "mapped_plus_strand_ref": ref_plus,
            "mapped_plus_strand_alt": alt_plus,
            "index_start": rel_start,
            "index_end": rel_end,
            "genomic_pos": f"{chrom}:{start}-{end}",
            "template_sequence": alt_template 
        },
        "results": results
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--alt", required=True)
    parser.add_argument("--strand", choices=["+", "-"], default="+", help="Strand of the input variant (+ or -)")
    parser.add_argument("--fasta", required=True)
    parser.add_argument("-p", "--padding", type=int, default=100)
    parser.add_argument("-k", "--top_k", type=int, default=5)
    parser.add_argument("--base_config", default="pcr/config/base_pcr.yaml")
    parser.add_argument("--system_config", default="pcr/config/system.yaml")
    
    args = parser.parse_args()
    
    res = design_qpcr_primers(
        args.chrom, args.start, args.end, args.ref, args.alt, args.strand,
        args.fasta, args.padding, args.top_k,
        args.base_config, args.system_config
    )
    print(json.dumps(res, indent=2))