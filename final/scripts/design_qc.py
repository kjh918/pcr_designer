#!/usr/bin/env python3
import sys
import os
import json
import argparse
from typing import Dict, Any, List, Optional

try:
    import primer3
except ImportError:
    primer3 = None

# 상위 경로 추가 (pcr 패키지 인식용)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import load_pipeline_config
from pcr.components.primer import Primer, Probe
from pcr.components.amplicon import Amplicon
from pcr.designers.base.qc import BaseQCExecutor
from pcr.qc.blast import BlastSpecificityChecker 


def reverse_complement(seq: str) -> str:
    """염기서열의 역상보 서열을 반환합니다."""
    return seq.translate(str.maketrans('ATGCatgcNn', 'TACGtacgNn'))[::-1]


def evaluate_qc_pipeline(
    project_name: str,
    fwd_seq: str,
    rev_seq: str,
    probe_seq: str = "",
    template_seq: str = "",
    genome: str = "hg38",
    qc_overrides: Optional[Dict[str, Any]] = None,
    base_yaml: str = "pcr/config/base_pcr.yaml", 
    system_yaml: str = "pcr/config/system.yaml"
) -> Dict[str, Any]:
    """
    웹(UI)에서 받은 서열 기반으로 열역학(BaseQC) 및 특이성(BlastQC) 검증을 통합 수행합니다.
    """
    assay_type = "qpcr"
    
    user_overrides = {}
    if qc_overrides:
        user_overrides["qc_criteria"] = qc_overrides

    config = load_pipeline_config(base_yaml, system_yaml, assay_type, user_overrides=user_overrides)
    
    amp_size = 0
    alignment_lines = []
    fwd_idx, rev_idx, prb_idx = 0, 0, 0
    
    if template_seq:
        fwd_idx = template_seq.find(fwd_seq)
        rev_rc = reverse_complement(rev_seq)
        rev_idx = template_seq.find(rev_rc)
        
        if fwd_idx != -1 and rev_idx != -1 and rev_idx >= fwd_idx:
            amp_size = rev_idx + len(rev_rc) - fwd_idx
            alignment_lines.append("[Input Template Alignment]")
            alignment_lines.append(f"TEMPLATE : {template_seq}")
            
            fwd_pad = " " * fwd_idx
            fwd_match = " " * fwd_idx + "|" * len(fwd_seq)
            alignment_lines.append(f"         : {fwd_match}")
            alignment_lines.append(f"FORWARD  : {fwd_pad}{fwd_seq}")
            
            rev_pad = " " * rev_idx
            rev_match = " " * rev_idx + "|" * len(rev_rc)
            alignment_lines.append(f"         : {rev_match}")
            alignment_lines.append(f"REV_RC   : {rev_pad}{rev_rc}")
            
            if probe_seq:
                prb_idx = template_seq.find(probe_seq)
                if prb_idx != -1:
                    prb_pad = " " * prb_idx
                    prb_match = " " * prb_idx + "|" * len(probe_seq)
                    alignment_lines.append(f"         : {prb_match}")
                    alignment_lines.append(f"PROBE    : {prb_pad}{probe_seq}")
                else:
                    prb_rc = reverse_complement(probe_seq)
                    prb_idx = template_seq.find(prb_rc)
                    if prb_idx != -1:
                        prb_pad = " " * prb_idx
                        prb_match = " " * prb_idx + "|" * len(prb_rc)
                        alignment_lines.append(f"         : {prb_match}")
                        alignment_lines.append(f"PROBE_RC : {prb_pad}{prb_rc}")

    fwd_primer = Primer.make_primer(
        sequence=fwd_seq,
        role="FORWARD",
        start_index=max(0, fwd_idx),
        end_index=max(0, fwd_idx) + len(fwd_seq)
    )
    
    rev_primer = Primer.make_primer(
        sequence=rev_seq,
        role="REVERSE",
        start_index=max(0, rev_idx),
        end_index=max(0, rev_idx) + len(rev_seq)
    )
    
    probe_obj = Probe.make_primer(
        sequence=probe_seq,
        role="INTERNAL",
        start_index=max(0, prb_idx),
        end_index=max(0, prb_idx) + len(probe_seq)
    ) if probe_seq else None

    virtual_template = template_seq if template_seq else fwd_seq + ("N" * 20) + reverse_complement(rev_seq)
    
    amp = Amplicon(
        id=f"{project_name}_1",
        set_id=project_name,
        forward=fwd_primer,
        reverse=rev_primer,
        probe=probe_obj,
        template_sequence=virtual_template
    )

    amp.product_size = amp_size if amp_size > 0 else len(virtual_template)

    blast_details = None
    blast_qc_pass = True
    blast_qc_log = ""
    
    if genome.lower() != "none":
        blast_checker = BlastSpecificityChecker(config)
        blast_checker.run([amp])
        blast_details = getattr(amp, "blast_stats", None)
        blast_qc_pass = getattr(amp, "is_qc_pass", True)
        blast_qc_log = getattr(amp, "qc_log", "")
        
        if blast_details:
            if blast_details.get("target_signals"):
                t_sig = blast_details["target_signals"][0]
                amp.product_size = t_sig["product_size"]
                amp.genomic_pos = t_sig["location"]

    try:
        qc_executor = BaseQCExecutor(config)
        evaluated_amps = qc_executor.execute([amp])
        amp = evaluated_amps[0]
    except Exception as e:
        amp.is_qc_pass = False
        amp.qc_log = f"BaseQC Error: {str(e)}"
        
    base_qc_pass = getattr(amp, "is_qc_pass", True)
    base_qc_log = getattr(amp, "qc_log", "")
    
    max_tm_diff = qc_overrides.get("oligo", {}).get("max_tm_diff", 3.0) if qc_overrides else 3.0
    tm_diff = abs(amp.forward.tm - amp.reverse.tm)
    if tm_diff > max_tm_diff:
        base_qc_pass = False
        base_qc_log += f" [Tm Diff {round(tm_diff,1)} > {max_tm_diff}]"

    amp.is_qc_pass = base_qc_pass and blast_qc_pass
    combined_logs = []
    if not base_qc_pass and base_qc_log: combined_logs.append(f"Thermo: {base_qc_log.strip()}")
    if not blast_qc_pass and blast_qc_log: combined_logs.append(f"BLAST: {blast_qc_log.strip()}")
    amp.qc_log = " | ".join(combined_logs) if combined_logs else "Passed: Thermally stable and Specific."

    # -------------------------------------------------------------------------
    # 🔥 [사용자 제안 로직 적용] 모든 앰플리콘을 독립적인 배열 객체로 분류
    # -------------------------------------------------------------------------
    results_list = []
    template_aln_str = "\n".join(alignment_lines)
    
    base_item = {
        "forward_primer": amp.forward.sequence,
        "reverse_primer": amp.reverse.sequence,
        "probe": amp.probe.sequence if amp.probe else "-",
        "tm_f": round(amp.forward.tm, 2),
        "tm_r": round(amp.reverse.tm, 2),
        "tm_p": round(amp.probe.tm, 2) if amp.probe else 0.0,
        "gc_f": round(getattr(amp.forward, "gc_percent", getattr(amp.forward, "gc", 0.0)), 2),
        "gc_r": round(getattr(amp.reverse, "gc_percent", getattr(amp.reverse, "gc", 0.0)), 2),
        "gc_p": round(getattr(amp.probe, "gc_percent", getattr(amp.probe, "gc", 0.0)), 2) if amp.probe else 0.0,
        "blast_stats": blast_details 
    }

    if blast_details and blast_details.get("total_signal_count", 0) > 0:
        rank = 1
        
        # 1. Target Signals (정상)
        for sig in blast_details.get("target_signals", []):
            item = base_item.copy()
            item["rank"] = rank
            item["id"] = f"{amp.id}_Target_{rank}"
            item["amplicon_size"] = sig.get("product_size", "N/A")
            item["genomic_pos"] = sig.get("location", "Unknown")
            
            # 메인 타겟일 경우 입력한 템플릿 매칭 뷰를 위에 덧붙임
            text_block = template_aln_str + ("\n\n" if template_aln_str else "") + sig.get("unified_text_block", "")
            item["alignment_text_block"] = text_block
            
            item["qc_info"] = {
                "is_pass": amp.is_qc_pass,
                "fail_reason": amp.qc_log if not amp.is_qc_pass else "PASS: Specific and stable target."
            }
            results_list.append(item)
            rank += 1
            
        # 2. Off-Target Signals (Probe 결합 포함된 비특이 앰플리콘)
        for sig in blast_details.get("off_target_signals", []):
            item = base_item.copy()
            item["rank"] = rank
            item["id"] = f"{amp.id}_OffTarget_{rank}"
            item["amplicon_size"] = sig.get("product_size", "N/A")
            item["genomic_pos"] = sig.get("location", "Unknown")
            item["alignment_text_block"] = sig.get("unified_text_block", "")
            item["qc_info"] = {
                "is_pass": False,
                "fail_reason": "FAIL: Off-Target Amplicon (Probe binds here!)"
            }
            results_list.append(item)
            rank += 1
            
        # 3. Amplification Only Signals (Probe 없는 단순 증폭 노이즈)
        for sig in blast_details.get("amplification_only", []):
            item = base_item.copy()
            item["rank"] = rank
            item["id"] = f"{amp.id}_Noise_{rank}"
            item["amplicon_size"] = sig.get("product_size", "N/A")
            item["genomic_pos"] = sig.get("location", "Unknown")
            item["alignment_text_block"] = sig.get("unified_text_block", "")
            item["qc_info"] = {
                "is_pass": False,
                "fail_reason": "WARNING: Amplification only (No probe binding)"
            }
            results_list.append(item)
            rank += 1

    else:
        # BLAST를 스킵했거나 Hit가 아예 없는 경우
        item = base_item.copy()
        item["rank"] = 1
        item["id"] = amp.id
        item["amplicon_size"] = amp.product_size if amp.product_size != len(virtual_template) else "N/A"
        item["genomic_pos"] = getattr(amp, "genomic_pos", "Unknown")
        item["alignment_text_block"] = template_aln_str
        
        fail_msg = amp.qc_log
        if blast_details and blast_details.get("total_signal_count", 0) == 0:
            fail_msg = "FAIL: No target found in reference genome."
            item["qc_info"] = {"is_pass": False, "fail_reason": fail_msg}
        else:
            item["qc_info"] = {
                "is_pass": getattr(amp, "is_qc_pass", True),
                "fail_reason": fail_msg if not getattr(amp, "is_qc_pass", True) else ""
            }
            
        results_list.append(item)

    # 모든 리스트를 통째로 전달
    return {
        "status": "success",
        "single_result": {"total_count": len(results_list)},
        "single_total_amplicons": results_list,
        "single_filtered_amplicons": results_list
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run QC Pipeline for Provided Sequences")
    parser.add_argument("--name", type=str, default="QC_Project", help="Project Name")
    parser.add_argument("--fwd", type=str, required=True, help="Forward Primer Sequence")
    parser.add_argument("--rev", type=str, required=True, help="Reverse Primer Sequence")
    parser.add_argument("--probe", type=str, default="", help="Probe Sequence (Optional)")
    parser.add_argument("--template", type=str, default="", help="Template Sequence (Optional)")
    parser.add_argument("--genome", type=str, default="hg38", help="Reference Genome for BLAST (hg38, mm10, none)")
    parser.add_argument("--base_config", type=str, default="pcr/config/base_pcr.yaml")
    parser.add_argument("--system_config", type=str, default="pcr/config/system.yaml")
    
    args = parser.parse_args()
    qc_overrides = {}

    try:
        res = evaluate_qc_pipeline(
            project_name=args.name, fwd_seq=args.fwd, rev_seq=args.rev, probe_seq=args.probe,
            template_seq=args.template, genome=args.genome, qc_overrides=qc_overrides,
            base_yaml=args.base_config, system_yaml=args.system_config
        )
        print(json.dumps(res, indent=2))
        
    except Exception as e:
        print(json.dumps({"status": "fail", "reason": str(e)}, indent=2))