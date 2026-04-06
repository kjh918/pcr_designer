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

from pcr.designers.base.schema import BaseDesignOutput
from pcr.designers.base.designer import BasePrimerDesigner
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

    # 1. 설정 로드
    config = load_pipeline_config(base_yaml, system_yaml, assay_type, user_overrides=user_overrides)

    # 🔥 [Error Fix] PCRParams 모델에 존재하지 않는 reference_name 필드 할당 시도를 제거합니다.
    # Pydantic strict mode로 인해 정의되지 않은 필드 할당 시 ValueError가 발생합니다.
    # 대신 하단의 BLAST 경로 설정을 통해 실행 여부를 제어합니다.

    # 2. BLAST 경로 설정 및 생략 제어

    if hasattr(config, "system") and hasattr(config.system, "paths"):
        if genome.lower() == "none":
            # BLAST를 돌리지 않으려면 DB 경로를 None으로 설정합니다.
            config.system.paths.blast_db_path = None
        else:
            # 유전체가 지정된 경우 해당 DB 경로 연결
            old_db = getattr(config.system.paths, "blast_db_path", "")
            db_dir = os.path.dirname(old_db) if old_db else f"db/{genome}"
            config.system.paths.blast_db_path = os.path.join(db_dir, genome)
    
    # 3. 서열 처리 및 정렬 시각화 데이터 생성 (기존 로직 유지)
    fwd_seq = fwd_seq.upper()
    rev_seq = rev_seq.upper()
    probe_seq = probe_seq.upper() if probe_seq else ""
    template_seq = template_seq.upper() if template_seq else ""

    amp_size = 0
    alignment_lines = []
    fwd_idx, rev_idx, prb_idx = -1, -1, -1

    display_template = template_seq
    
    if template_seq:
        fwd_idx = template_seq.find(fwd_seq)
        rev_rc = reverse_complement(rev_seq)
        rev_idx = template_seq.find(rev_rc)
        
        if '[' in template_seq and ']' in template_seq:
            clean_tpl, ref_tpl, target_indices = BasePrimerDesigner.parse_sequence_with_brackets(template_seq)
            work_template = clean_tpl
        
        fwd_idx = work_template.find(fwd_seq)
        rev_rc = reverse_complement(rev_seq)
        rev_idx = work_template.find(rev_rc)
        
        if fwd_idx != -1 and rev_idx != -1 and rev_idx >= fwd_idx:
            amp_size = rev_idx + len(rev_rc) - fwd_idx
            alignment_lines.append("[Input Template Alignment]")
            alignment_lines.append(f"TEMPLATE: {ref_tpl}")
            
            # 정렬 가이드 라인은 가변적인 대괄호 길이를 고려하기 위해 clean_seq(work_template) 기준으로 계산하되
            # 시각화는 텍스트 블록의 특성상 clean 서열을 보여주는 것이 더 정확하므로 라벨을 업데이트합니다.
            if display_template != work_template:
                alignment_lines.append(f"(ALT)   : {work_template}")

            alignment_lines.append(f"        : " + " " * fwd_idx + "|" * len(fwd_seq))
            alignment_lines.append(f"FORWARD : " + " " * fwd_idx + fwd_seq)
            alignment_lines.append(f"        : " + " " * rev_idx + "|" * len(rev_rc))
            alignment_lines.append(f"REV_RC  : " + " " * rev_idx + rev_rc)
            print(work_template)
            print(probe_seq)

            if probe_seq:
                prb_idx = work_template.find(probe_seq)
                if prb_idx == -1:
                    prb_rc = reverse_complement(probe_seq)
                    prb_idx = work_template.find(prb_rc)
                if prb_idx != -1:
                    alignment_lines.append(f"        : " + " " * prb_idx + "|" * len(probe_seq))
                    alignment_lines.append(f"PROBE   : " + " " * prb_idx + (probe_seq if work_template.find(probe_seq) != -1 else prb_rc))

    # 4. 앰플리콘 객체 생성
    fwd_primer = Primer.make_primer(sequence=fwd_seq, role="FORWARD", start_index=max(0, fwd_idx), end_index=max(0, fwd_idx) + len(fwd_seq))
    rev_primer = Primer.make_primer(sequence=rev_seq, role="REVERSE", start_index=max(0, rev_idx), end_index=max(0, rev_idx) + len(rev_seq))
    probe_obj = Probe.make_primer(sequence=probe_seq, role="INTERNAL", start_index=max(0, prb_idx), end_index=max(0, prb_idx) + len(probe_seq)) if probe_seq else None

    virtual_template = template_seq if template_seq else fwd_seq + ("N" * 20) + reverse_complement(rev_seq)
    
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

    try:
        qc_executor = QPCRQCExecutor(config,reference_name=genome.lower())
        evaluated_amps = qc_executor.execute([amp])
        
        output_schema = BaseDesignOutput(
            status="success",
            amplicons=evaluated_amps
        )
        
        return output_schema.to_frontend_dict()
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"status": "fail", "reason": f"QC 파이프라인 오류: {str(e)}"}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", type=str, default="QC_Project")
    parser.add_argument("--fwd", type=str, required=True)
    parser.add_argument("--rev", type=str, required=True)
    parser.add_argument("--probe", type=str, default="")
    parser.add_argument("--template", type=str, default="")
    parser.add_argument("--genome", type=str, default="hg38")
    parser.add_argument("--base_config", type=str, default="pcr/config/base_pcr.yaml")
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