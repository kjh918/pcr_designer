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
from pcr.designers.base.designer import BasePrimerDesigner
from pcr.designers.base.schema import BaseDesignOutput

# 🔥 최신 API 스키마 임포트
from pcr.designers.qpcr.schema import QPCRDesignInput

def design_qpcr_pipeline(
    req: QPCRDesignInput,
    designer_yaml: str = "pcr/designers/qpcr/config.yaml", 
    system_yaml: str = "pcr/config/system.yaml"
) -> Dict[str, Any]:
    """
    최신 Pydantic 스키마(req)를 통째로 받아 코어 엔진을 구동하고,
    규격화된 BaseDesignOutput을 통해 결과를 반환하는 파이프라인입니다.
    """
    assay_type = "qpcr"

    # 🔥 타겟 자동 파싱: 프론트에서 대괄호([ ])로 감싸 보낸 서열을 자동으로 분해
    clean_seq, target_indices = BasePrimerDesigner.parse_sequence_with_brackets(req.sequence)
    
    if not target_indices:
        return {"status": "fail", "reason": "Target not found. Please wrap your target with brackets, e.g., ATGC[A/G]ATGC"}

    rel_start = min(target_indices) - 1
    rel_end = max(target_indices)

    # 1. 프론트엔드에서 넘어온 딕셔너리 파라미터를 오버라이드용으로 조립
    user_overrides = {
        "pcr_params": req.pcr_params,
        "qc_criteria": req.qc_criteria
    }

    # config loader 실행
    config = load_pipeline_config(designer_yaml, system_yaml, assay_type, user_overrides=user_overrides)

    # 2. [In-Memory Override] BLAST 시스템 경로 동적 할당
    if hasattr(config, "system") and hasattr(config.system, "paths"):
        if req.reference_genome.lower() == "none":
            config.system.paths.blast_db_path = None
        else:
            old_db = getattr(config.system.paths, "blast_db_path", "")
            db_dir = os.path.dirname(old_db) if old_db else f"db/{req.reference_genome}"
            config.system.paths.blast_db_path = os.path.join(db_dir, req.reference_genome)

    # 3. Factory 실행
    factory = PCRFactory(config)

    output = factory.run(
        assay_type=assay_type,
        name=req.design_name,
        target_start=rel_start,
        target_end=rel_end,
        reference_name=req.reference_genome,
        template_sequence=clean_seq,
        top_k=req.top_k,
        run_qc=True,
        # 매뉴얼 서열 입력이므로 genomic start와 strand를 기본값으로 고정하여 안정성 확보
        template_genomic_start=0,
        target_strand="+"
    )

    if output.status != "success":
        return {
            "status": "fail", 
            "reason": getattr(output, "error_msg", "Design failed"),
            "log_messages": getattr(output, "log_messages", "")
        }

    # 4. 공통 스키마(BaseDesignOutput)를 활용한 압도적으로 깔끔한 결과 자동 포맷팅
    # 기존 코드에 있던 수십 줄의 딕셔너리 매핑과 get_position_info()가 이 한 줄로 대체됩니다!
    output_schema = BaseDesignOutput(
        status="success",
        amplicons=output.amplicons
    )
    
    result_dict = output_schema.to_frontend_dict()
    
    # 라우터에서 Export 메타데이터로 사용할 수 있도록 원본 파라미터 백업
    result_dict["_applied_pcr_params"] = req.pcr_params
    result_dict["_applied_qc_params"] = req.qc_criteria
    
    return result_dict


if __name__ == "__main__":
    # 단독 실행(CLI) 테스트용 안전 장치
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", type=str, default="CLI_qPCR_Project")
    parser.add_argument("--seq", type=str, required=True, help="Sequence with brackets, e.g., ATGC[A]ATGC")
    parser.add_argument("--genome", type=str, default="none")
    args = parser.parse_args()

    # CLI에서도 QPCRDesignApiInput 스키마를 목업하여 전달
    mock_req = QPCRDesignApiInput(
        design_name=args.name,
        sequence=args.seq,
        reference_genome=args.genome
    )

    try:
        res = design_qpcr_pipeline(mock_req)
        print(json.dumps(res, indent=2)) 
    except Exception as e:
        print(json.dumps({"status": "fail", "reason": str(e)}, indent=2))