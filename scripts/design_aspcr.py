#!/usr/bin/env python3
import sys
import os
import json
import argparse
from typing import Dict, Any, Optional

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import load_pipeline_config
from pcr.factory import PCRFactory
from pcr.designers.as_pcr.schema import ASPCRDesignOutput

def _generate_aspcr_templates(
    left_seq: str, right_seq: str, ref: str, alt: str,
    fixed_prime: str = "forward", 
    mismatch_pos: int = 3,
    mismatch_intensity: str = "strong"
) -> Dict[str, str]:
    def _apply_as_pcr_logic(target_base: str, intensity: str) -> str:
        intensity=intensity.lower()
        mismatch_map = {
            "A": "G" if intensity == "strong" else "C",
            "G": "A" if intensity == "strong" else "T",
            "C": "T" if intensity == "strong" else "A",
            "T": "C" if intensity == "strong" else "G"
        }
        return mismatch_map.get(target_base.upper(), "N")

    wt_mm_left, alt_mm_left = left_seq, left_seq
    wt_mm_right, alt_mm_right = right_seq, right_seq

    if mismatch_pos and mismatch_pos > 1:
        if fixed_prime == "forward":
            mm_index = -(mismatch_pos - 1)
            if len(left_seq) >= abs(mm_index):
                original_base = left_seq[mm_index]
                new_base = _apply_as_pcr_logic(original_base, mismatch_intensity)
                wt_mm_left = left_seq[:mm_index] + new_base + left_seq[mm_index+1:]
                alt_mm_left = wt_mm_left
        else: # reverse
            mm_index = mismatch_pos - 2
            if len(right_seq) > mm_index:
                original_base = right_seq[mm_index]
                new_base = _apply_as_pcr_logic(original_base, mismatch_intensity)
                wt_mm_right = right_seq[:mm_index] + new_base + right_seq[mm_index+1:]
                alt_mm_right = wt_mm_right

    return {
        "wt": left_seq + ref + right_seq,
        "alt": left_seq + alt + right_seq,
        "wt_mm": wt_mm_left + ref + wt_mm_right,
        "alt_mm": alt_mm_left + alt + alt_mm_right
    }

def design_aspcr_primers(
    design_name: str,
    sequence: str,
    genome: str = "hg38",
    top_k: int = 5,
    base_yaml: str = "pcr/config/base_pcr.yaml", 
    system_yaml: str = "pcr/config/system.yaml",
    pcr_overrides: Optional[Dict[str, Any]] = None,
    qc_overrides: Optional[Dict[str, Any]] = None,
    fixed_prime: str = "forward", 
    mismatch_pos: int = 3,
    mismatch_intensity: str = "strong"
) -> Dict[str, Any]:
    
    assay_type = "aspcr"

    start_idx = sequence.find('[')
    end_idx = sequence.find(']')
    if start_idx == -1 or end_idx == -1:
        return {"status": "fail", "reason": "Target brackets [REF,ALT] not found."}

    left_seq = sequence[:start_idx].replace(" ", "").upper()
    right_seq = sequence[end_idx+1:].replace(" ", "").upper()
    alleles = sequence[start_idx+1:end_idx].replace(" ", "").upper()
    
    if ',' in alleles:
        ref, alt = alleles.split(',')
    elif '/' in alleles:
        ref, alt = alleles.split('/')
    else:
        ref, alt = alleles, alleles

    print(fixed_prime)
    print(mismatch_pos)
    print(mismatch_intensity)
    
    templates = _generate_aspcr_templates(
        left_seq, right_seq, ref, alt, 
        fixed_prime=fixed_prime,
        mismatch_pos=mismatch_pos,
        mismatch_intensity=mismatch_intensity
    )

    user_overrides = {}
    if qc_overrides: user_overrides["qc_criteria"] = qc_overrides
    if pcr_overrides: user_overrides["pcr_params"] = pcr_overrides

    config = load_pipeline_config(base_yaml, system_yaml, assay_type, user_overrides=user_overrides)
    
    if hasattr(config, "system") and hasattr(config.system, "paths"):
        if genome.lower() == "none":
            config.system.paths.blast_db_path = None
        else:
            old_db = getattr(config.system.paths, "blast_db_path", "")
            db_dir = os.path.dirname(old_db) if old_db else f"db/{genome}"
            config.system.paths.blast_db_path = os.path.join(db_dir, genome)

    factory = PCRFactory(config)
    
    output = factory.run(
        assay_type=assay_type,
        name=design_name,
        reference_name=genome,
        template_sequence=sequence,
        top_k=top_k,
        run_qc=True,
        overrides={"PRIMER_NUM_RETURN": 30},
        templates=templates,
        fixed_prime=fixed_prime,
        mismatch_pos=mismatch_pos,
        mismatch_intensity=mismatch_intensity
    )

    if output.status != "success":
        fail_reason = output.error_msg or "Failed to design AS-PCR primers."
        logs = getattr(output, "log_messages", [])
        return {
            "status": "fail", 
            "reason": fail_reason,
            "log_messages": "; ".join(logs) if isinstance(logs, list) else str(logs)
        }

    final_output = ASPCRDesignOutput(
        status=output.status,
        amplicons=output.amplicons,
        error_msg=output.error_msg,
        log_messages=output.log_messages,
        metadata={
            "assay": "AS-PCR", 
            "design_name": design_name,
            "genome_build": genome,
            "fixed_prime": fixed_prime,
            "mismatch_pos": mismatch_pos,
            "mismatch_intensity": mismatch_intensity
        }
    )
    
    # qc.py 계층과 동일한 포맷으로 래핑하여 리턴
    frontend_dict = final_output.to_frontend_dict()
    results_list = frontend_dict.get("results", [])
    
    return {
        "status": "success",
        "metadata": {
            "project_name": design_name,
            "reference_genome": genome,
            "assay": "AS-PCR"
        },
        "summary": frontend_dict.get("summary", {}),
        "inputs": {
            "pcr_params": pcr_overrides or {},
            "qc_criteria": qc_overrides or {}
        },
        "results": results_list
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Design AS-PCR primers via Factory")
    parser.add_argument("--name", default="AS-PCR_Test")
    parser.add_argument("--seq", required=True, help="Sequence with brackets, e.g., ATGC[A,G]ATGC")
    parser.add_argument("--genome", default="none")
    parser.add_argument("-k", "--top_k", type=int, default=5)
    parser.add_argument("--fixed_prime", choices=["forward", "reverse"], default="forward")
    parser.add_argument("-m", "--mismatch_pos", type=int, default=3)
    parser.add_argument("--intensity", choices=["strong", "weak"], default="strong")
    parser.add_argument("--base_config", default="pcr/config/base_pcr.yaml")
    parser.add_argument("--system_config", default="pcr/config/system.yaml")
    
    args = parser.parse_args()
    
    res = design_aspcr_primers(
        design_name=args.name,
        sequence=args.seq,
        genome=args.genome,
        top_k=args.top_k,
        base_yaml=args.base_config,
        system_yaml=args.system_config,
        fixed_prime=args.fixed_prime,
        mismatch_pos=args.mismatch_pos,       
        mismatch_intensity=args.intensity     
    )
    
    print(json.dumps(res, indent=2))