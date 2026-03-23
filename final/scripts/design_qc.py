#!/usr/bin/env python3
import sys
import os
import json
import argparse
from typing import Dict, Any, Optional

try:
    import primer3
except ImportError:
    primer3 = None

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import load_pipeline_config
from pcr.components.primer import Primer, Probe
from pcr.components.amplicon import Amplicon

# 🔥 복잡한 수동 BLAST 로직은 모두 제거하고, 
# 똑똑해진 QPCRQCExecutor와 스키마만 가져옵니다!
from pcr.designers.base.schema import BaseDesignOutput
from pcr.designers.qpcr.qc import QPCRQCExecutor


def reverse_complement(seq: str) -> str:
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
    
    assay_type = "qpcr"
    
    user_overrides = {}
    if qc_overrides:
        user_overrides["qc_criteria"] = qc_overrides

    config = load_pipeline_config(base_yaml, system_yaml, assay_type, user_overrides=user_overrides)
    
    # [In-Memory Override] BLAST 시스템 경로 동적 할당
    if hasattr(config, "system") and hasattr(config.system, "paths"):
        if genome.lower() == "none":
            config.system.paths.blast_db_path = None
        else:
            old_db = getattr(config.system.paths, "blast_db_path", "")
            db_dir = os.path.dirname(old_db) if old_db else f"db/{genome}"
            config.system.paths.blast_db_path = os.path.join(db_dir, genome)

    fwd_seq = fwd_seq.upper()
    rev_seq = rev_seq.upper()
    probe_seq = probe_seq.upper() if probe_seq else ""
    template_seq = template_seq.upper() if template_seq else ""

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
            alignment_lines.append(f"         : " + " " * fwd_idx + "|" * len(fwd_seq))
            alignment_lines.append(f"FORWARD  : " + " " * fwd_idx + fwd_seq)
            alignment_lines.append(f"         : " + " " * rev_idx + "|" * len(rev_rc))
            alignment_lines.append(f"REV_RC   : " + " " * rev_idx + rev_rc)
            
            if probe_seq:
                prb_idx = template_seq.find(probe_seq)
                if prb_idx != -1:
                    alignment_lines.append(f"         : " + " " * prb_idx + "|" * len(probe_seq))
                    alignment_lines.append(f"PROBE    : " + " " * prb_idx + probe_seq)
                else:
                    prb_rc = reverse_complement(probe_seq)
                    prb_idx = template_seq.find(prb_rc)
                    if prb_idx != -1:
                        alignment_lines.append(f"         : " + " " * prb_idx + "|" * len(prb_rc))
                        alignment_lines.append(f"PROBE_RC : " + " " * prb_idx + prb_rc)

    fwd_primer = Primer.make_primer(sequence=fwd_seq, role="FORWARD", start_index=max(0, fwd_idx), end_index=max(0, fwd_idx) + len(fwd_seq))
    rev_primer = Primer.make_primer(sequence=rev_seq, role="REVERSE", start_index=max(0, rev_idx), end_index=max(0, rev_idx) + len(rev_seq))
    probe_obj = Probe.make_primer(sequence=probe_seq, role="INTERNAL", start_index=max(0, prb_idx), end_index=max(0, prb_idx) + len(probe_seq)) if probe_seq else None

    virtual_template = template_seq if template_seq else fwd_seq + ("N" * 20) + reverse_complement(rev_seq)
    
    # 💡 1. 스크립트의 역할: 사용자가 준 서열로 "단 1개의 Amplicon"만 조립합니다.
    amp = Amplicon(
        id=f"{project_name}_1", set_id=project_name,
        forward=fwd_primer, reverse=rev_primer, probe=probe_obj,
        template_sequence=virtual_template
    )
    amp.product_size = amp_size if amp_size > 0 else len(virtual_template)
    amp.alignment_visual = alignment_lines

    try:
        amp.gc_percent = (virtual_template.count('G') + virtual_template.count('C')) / len(virtual_template) * 100 if virtual_template else 0.0
        amp.tm = primer3.calc_tm(virtual_template) if primer3 and virtual_template else 0.0
    except Exception:
        pass


    # =====================================================================
    # 💡 2. 아키텍처의 꽃: 파이프라인(Executor) 호출 및 Schema 패키징
    # =====================================================================
    try:
        # 단 1줄! QPCRQCExecutor가 내부적으로 BLAST를 돌리고, 
        # Hit가 5개라면 5개로 증식(Explode)시킨 앰플리콘 리스트를 반환합니다.
        qc_executor = QPCRQCExecutor(config)
        evaluated_amps = qc_executor.execute([amp])
        
        # 증식되어 검증까지 완료된 리스트를 스키마에 던져 넣으면 프론트용 JSON이 튀어나옵니다.
        output_schema = BaseDesignOutput(
            status="success",
            amplicons=evaluated_amps
        )
        
        return output_schema.to_frontend_dict()
        
    except Exception as e:
        return {"status": "fail", "reason": f"QC 파이프라인 오류: {str(e)}"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", type=str, default="QC_Project")
    parser.add_argument("--fwd", type=str, required=True)
    parser.add_argument("--rev", type=str, required=True)
    parser.add_argument("--probe", type=str, default="")
    parser.add_argument("--template", type=str, default="")
    parser.add_argument("--genome", type=str, default="hg38")
    parser.add_argument("--base_config", type=str, default="pcr/designers/qpcr/config.yaml")
    parser.add_argument("--system_config", type=str, default="pcr/config/system.yaml")
    args = parser.parse_args()

    try:
        res = evaluate_qc_pipeline(
            project_name=args.name, fwd_seq=args.fwd, rev_seq=args.rev, probe_seq=args.probe,
            template_seq=args.template, genome=args.genome, qc_overrides={},
            base_yaml=args.base_config, system_yaml=args.system_config
        )
        print(json.dumps(res, indent=2)) 
    except Exception as e:
        print(json.dumps({"status": "fail", "reason": str(e)}, indent=2))