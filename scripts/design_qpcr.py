#!/usr/bin/env python3
import sys
import os
import json
import argparse
from typing import Dict, Any

# 상위 경로 추가 (pcr 패키지 인식용)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import load_pipeline_config
from pcr.factory import PCRFactory
from pcr.designers.qpcr.schema import QPCRDesignOutput

def design_qpcr_pipeline(
    design_name: str,
    sequence: str,
    genome: str = "hg38",
    top_k: int = 10,
    pcr_overrides: Dict[str, Any] = None,
    qc_overrides: Dict[str, Any] = None,
    designer_yaml: str = "pcr/designers/qpcr/config.yaml", 
    system_yaml: str = "pcr/config/system.yaml"
) -> Dict[str, Any]:
    """
    라우터가 풀어준 파라미터(kwargs)를 그대로 받아,
    PCRFactory 하나에 모든 과정(Design -> QC -> Ranking)을 위임(Delegation)합니다.
    (대괄호 파싱 역시 Factory가 내부적으로 자동 처리합니다.)
    """
    assay_type = "qpcr"

    # 1. Config 로드 및 오버라이드 
    user_overrides = {}
    if pcr_overrides:
        user_overrides["pcr_params"] = pcr_overrides
    if qc_overrides:
        user_overrides["qc_criteria"] = qc_overrides

    config = load_pipeline_config(designer_yaml, system_yaml, assay_type, user_overrides=user_overrides)
    
    if genome.lower() == "none" and hasattr(config, "system") and hasattr(config.system, "paths"):
        config.system.paths.blast_db_path = None
    
    reference_name = genome.lower()
    # 3. 🔥 PCR Factory 구동 (디자인부터 검증, 랭킹, 서열 파싱까지 Factory가 100% 처리합니다)
    factory = PCRFactory(config)
    
    try:
        output = factory.run(
            assay_type=assay_type,
            name=design_name,
            reference_name=genome,
            template_sequence=sequence,  # 🔥 괄호가 포함된 원본 서열을 그대로 던집니다!
            top_k=top_k,
            run_qc=True,
            template_genomic_start=0,
            target_strand="+"
        )
    except Exception as e:
        # 괄호가 없거나 파싱 오류 시 팩토리가 뱉는 에러를 우아하게 잡아냅니다.
        return {
            "status": "fail",
            "reason": f"Factory Pipeline Error: {str(e)}",
            "log_messages": str(e)
        }

    # 4. 에러 처리 및 반환
    if output.status != "success":
        fail_reason = getattr(output, "error_msg", "Failed to design primers/probes.")
        logs = getattr(output, "log_messages", [])
        return {
            "status": "fail", 
            "reason": fail_reason,
            "log_messages": "; ".join(logs) if isinstance(logs, list) else str(logs)
        }

    # 5. 결과를 순수 코어 모델인 QPCRDesignOutput에 얹어서 프론트엔드로 변환
    output_schema = QPCRDesignOutput(
        status="success",
        amplicons=output.amplicons,
        log_messages=output.log_messages
    )
    
    return output_schema.to_frontend_dict()


if __name__ == "__main__":
    # CLI 단독 실행용
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", type=str, default="CLI_qPCR_Project")
    parser.add_argument("--seq", type=str, required=True, help="Sequence with brackets, e.g., ATGC[A/G]ATGC")
    parser.add_argument("--genome", type=str, default="none")
    args = parser.parse_args()

    try:
        res = design_qpcr_pipeline(
            design_name=args.name,
            sequence=args.seq,
            genome=args.genome
        )
        print(json.dumps(res, indent=2)) 
    except Exception as e:
        print(json.dumps({"status": "fail", "reason": str(e)}, indent=2))