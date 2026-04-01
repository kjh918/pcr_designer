import os
import traceback
from datetime import datetime
from typing import Dict, Any

from fastapi import APIRouter, HTTPException

from pcr.designers.qpcr.schema import QPCRDesignInput
from scripts.design_qpcr import design_qpcr_pipeline 

router = APIRouter(
    prefix="/api/design",
    tags=["qpcr"]
)

@router.post("/qpcr")
async def design_qpcr_api(req: QPCRDesignInput): # 🔥 수정 2: 타입 힌트 변경
    print(f"\n🚀 [API] qPCR Design Request: {req.design_name}")
    
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__)) 
    APP_DIR = os.path.dirname(CURRENT_DIR) 
    ROOT_DIR = os.path.dirname(APP_DIR) 

    SYSTEM_YAML_PATH = os.path.join(ROOT_DIR, "pcr", "config", "system.yaml")
    DESIGNER_YAML_PATH = os.path.join(ROOT_DIR, "pcr", "designers", "base", "config.yaml")
    print(DESIGNER_YAML_PATH)
    try:
        # 🔥 수정 3: 수십 줄의 지저분한 수동 매핑 로직을 전부 지우고, 
        # schema.py에 만들어둔 우아한 Adapter 메서드를 바로 호출하여 딕셔너리를 뽑아냅니다.
        pcr_core_overrides = req.to_core_pcr_params()
        qc_core_overrides = req.to_core_qc_overrides()

        # 스크립트 실행
        raw_result = design_qpcr_pipeline(
            design_name=req.design_name,
            sequence=req.sequence,
            genome=req.reference_genome,
            top_k=req.top_k,
            pcr_overrides=pcr_core_overrides,
            qc_overrides=qc_core_overrides,
            designer_yaml=DESIGNER_YAML_PATH,
            system_yaml=SYSTEM_YAML_PATH
        )
        
        # 최종 출력 포맷 (qc.py와 완벽히 동일한 평탄화 계층)
        results_list = raw_result.get("results", [])
        
        final_output = {
            "status": raw_result.get("status", "error"),
            "metadata": {
                "project_name": req.design_name,
                "reference_genome": req.reference_genome,
                "timestamp": datetime.now().isoformat(timespec="seconds")
            },
            "summary": raw_result.get("summary", {}),
            "inputs": {
                "pcr_params": pcr_core_overrides,
                "qc_criteria": qc_core_overrides
            },
            "results": results_list
        }
        
        if final_output["status"] == "success":
            passed = final_output.get("summary", {}).get("passed_count", 0)
            total = final_output.get("summary", {}).get("total_count", 0)
            print(f"✅ qPCR Design Success. Passed {passed} / {total}")
            
        return final_output

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))