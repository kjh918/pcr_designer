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
    (템플릿이 없는 경우 BLAST를 통해 앰플리콘 사이즈와 위치를 역추적합니다.)
    """
    assay_type = "qpcr"
    
    # -------------------------------------------------------------------------
    # 1. Config 로드 및 BLAST 제어
    # -------------------------------------------------------------------------
    user_overrides = {}
    if qc_overrides:
        user_overrides["qc_criteria"] = qc_overrides

    config = load_pipeline_config(base_yaml, system_yaml, assay_type, user_overrides=user_overrides)
    
    # -------------------------------------------------------------------------
    # 2. 템플릿 매칭 (사전 위치 파악 - 템플릿이 주어진 경우에만)
    # -------------------------------------------------------------------------
    amp_size = 0
    alignment_lines = []
    fwd_idx, rev_idx, prb_idx = 0, 0, 0
    
    if template_seq:
        fwd_idx = template_seq.find(fwd_seq)
        rev_rc = reverse_complement(rev_seq)
        rev_idx = template_seq.find(rev_rc)
        
        if fwd_idx != -1 and rev_idx != -1 and rev_idx >= fwd_idx:
            amp_size = rev_idx + len(rev_rc) - fwd_idx
            alignment_lines.append("[Input Template Matching]")
            alignment_lines.append(f"TEMPLATE : {template_seq}")
            alignment_lines.append(f"FORWARD  : {' ' * fwd_idx}{fwd_seq}")
            alignment_lines.append(f"REVERSE_RC: {' ' * rev_idx}{rev_rc}")
            
            if probe_seq:
                prb_idx = template_seq.find(probe_seq)
                if prb_idx != -1:
                    alignment_lines.append(f"PROBE    : {' ' * prb_idx}{probe_seq}")
                else:
                    prb_rc = reverse_complement(probe_seq)
                    prb_idx = template_seq.find(prb_rc)
                    if prb_idx != -1:
                        alignment_lines.append(f"PROBE_RC : {' ' * prb_idx}{prb_rc}")

    # -------------------------------------------------------------------------
    # 3. 임시 객체 생성 (열역학 프로퍼티 자동 정의)
    # -------------------------------------------------------------------------
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
    )
    # 템플릿이 없으면 임시 가상 템플릿 생성
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
    #amp.reference_name = genome if genome.lower() != "none" else None

    # -------------------------------------------------------------------------
    # 4. BLAST 특이성 검사 (가장 먼저 실행하여 진짜 Amplicon Size 및 위치 획득)
    # -------------------------------------------------------------------------
    blast_details = None
    blast_qc_pass = True
    blast_qc_log = ""
    

    blast_checker = BlastSpecificityChecker(config)
    blast_checker.run([amp]) # 내부적으로 amp.blast_stats 및 is_qc_pass 업데이트됨
    blast_details = getattr(amp, "blast_stats", None)
    blast_qc_pass = getattr(amp, "is_qc_pass", True)
    blast_qc_log = getattr(amp, "qc_log", "")


    print(blast_details)
    print(blast_qc_log)
    
    # 🔥 BLAST에서 진짜 타겟(On-target)을 찾은 경우 Amplicon 객체 업데이트
    if blast_details and blast_details.get("target_signals"):
        t_sig = blast_details["target_signals"][0]
        
        # 템플릿 없이도 진짜 증폭 사이즈와 위치를 확보함
        amp.product_size = t_sig["product_size"]
        amp.genomic_pos = t_sig["location"]
        
        # 시각화 텍스트블록 처리
        if not template_seq:
            alignment_lines = ["[BLAST On-Target Alignment (Derived from Reference Genome)]"]
        else:
            alignment_lines.append("\n[BLAST On-Target Alignment]")
            
        alignment_lines.append(f"Location: {t_sig['location']} (Size: {t_sig['product_size']}bp)")
        
        for key in ["fwd", "rev", "probe"]:
            if t_sig.get(key):
                blk = t_sig[key]
                alignment_lines.append(f"\n--- {blk['label']} ({blk['identity']}) ---")
                alignment_lines.append(f"Q: {blk['query']}")
                alignment_lines.append(f"   {blk['match']}")
                alignment_lines.append(f"S: {blk['subject']}")
    try:
        qc_executor = BaseQCExecutor(config)
        evaluated_amps = qc_executor.execute([amp])
        amp = evaluated_amps[0]
    except Exception as e:
        amp.is_qc_pass = False
        amp.qc_log = f"BaseQC Error: {str(e)}"
        
    base_qc_pass = getattr(amp, "is_qc_pass", True)
    base_qc_log = getattr(amp, "qc_log", "")
    
    # 수동 Tm Diff 검사 추가
    max_tm_diff = qc_overrides.get("oligo", {}).get("max_tm_diff", 3.0) if qc_overrides else 3.0
    tm_diff = abs(amp.forward.tm - amp.reverse.tm)
    if tm_diff > max_tm_diff:
        base_qc_pass = False
        base_qc_log += f" [Tm Diff {round(tm_diff,1)} > {max_tm_diff}]"

    # -------------------------------------------------------------------------
    # 6. 최종 판정 (BaseQC와 BlastQC 병합)
    # -------------------------------------------------------------------------
    amp.is_qc_pass = base_qc_pass and blast_qc_pass
    combined_logs = []
    if not base_qc_pass and base_qc_log: combined_logs.append(f"Thermo: {base_qc_log.strip()}")
    if not blast_qc_pass and blast_qc_log: combined_logs.append(f"BLAST: {blast_qc_log.strip()}")
    amp.qc_log = " | ".join(combined_logs) if combined_logs else "Passed: Thermally stable and Specific."

    # -------------------------------------------------------------------------
    # 7. 프론트엔드 호환 포맷팅 반환
    # -------------------------------------------------------------------------
    final_amp_size = amp.product_size if amp.product_size != len(virtual_template) else "N/A"

    formatted_item = {
        "rank": 1,
        "id": amp.id,
        "forward_primer": amp.forward.sequence,
        "reverse_primer": amp.reverse.sequence,
        "probe": amp.probe.sequence if amp.probe else "-",
        "tm_f": round(amp.forward.tm, 2),
        "tm_r": round(amp.reverse.tm, 2),
        "tm_p": round(amp.probe.tm, 2) if amp.probe else 0.0,
        "gc_f": round(getattr(amp.forward, "gc_percent", getattr(amp.forward, "gc", 0.0)), 2),
        "gc_r": round(getattr(amp.reverse, "gc_percent", getattr(amp.reverse, "gc", 0.0)), 2),
        "gc_p": round(getattr(amp.probe, "gc_percent", getattr(amp.probe, "gc", 0.0)), 2) if amp.probe else 0.0,
        "amplicon_size": final_amp_size,
        "genomic_pos": getattr(amp, "genomic_pos", "Unknown"),
        "alignment_text_block": "\n".join(alignment_lines),
        "qc_info": {
            "is_pass": getattr(amp, "is_qc_pass", True),
            "fail_reason": getattr(amp, "qc_log", "") if not getattr(amp, "is_qc_pass", True) else ""
        },
        "blast_stats": blast_details # 프론트에서 Off-target 개수 등을 시각화할 때 활용 가능
    }

    return {
        "status": "success",
        "single_result": {"total_count": 1},
        "single_total_amplicons": [formatted_item],
        "single_filtered_amplicons": [formatted_item] if getattr(amp, "is_qc_pass", True) else []
    }

# =====================================================================
# CLI 실행부 (테스트용)
# =====================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run QC Pipeline for Provided Sequences")
    parser.add_argument("--name", type=str, default="QC_Project", help="Project Name")
    parser.add_argument("--fwd", type=str, required=True, help="Forward Primer Sequence")
    parser.add_argument("--rev", type=str, required=True, help="Reverse Primer Sequence")
    parser.add_argument("--probe", type=str, default="", help="Probe Sequence (Optional)")
    parser.add_argument("--template", type=str, default="", help="Template Sequence (Optional)")
    parser.add_argument("--genome", type=str, default="none", help="Reference Genome for BLAST (hg38, mm10, none)")
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