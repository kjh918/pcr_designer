import os
import traceback
from datetime import datetime
from fastapi import APIRouter, HTTPException

from pcr.designers.as_pcr.schema import ASPCRDesignInput
from scripts.design_aspcr import design_aspcr_primers 

router = APIRouter(
    prefix="/api/design",
    tags=["aspcr"]
)

@router.post("/aspcr")
async def design_aspcr_api(req: ASPCRDesignInput):
    print(f"\n🚀 [API] ASPCR Design Request: {req.design_name}")
    
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__)) 
    APP_DIR = os.path.dirname(CURRENT_DIR) 
    ROOT_DIR = os.path.dirname(APP_DIR) 

    SYSTEM_YAML_PATH = os.path.join(ROOT_DIR, "pcr", "config", "system.yaml")
    ASPCR_CONFIG_PATH = os.path.join(ROOT_DIR, "pcr", "designers", "as_pcr", "config.yaml")
    try:
        # 🔥 Step 1: Pydantic 스키마를 통해 백엔드 코어용 Dict 자동 조립
        pcr_overrides_dict = req.to_core_pcr_params()
        qc_overrides_dict = req.to_core_qc_overrides()

        # 🔥 Step 2: 코어 스크립트 실행
        raw_result = design_aspcr_primers(
            design_name=req.design_name,
            sequence=req.sequence,
            genome=req.reference_genome,
            top_k=req.top_k,
            base_yaml=ASPCR_CONFIG_PATH,
            system_yaml=SYSTEM_YAML_PATH,
            pcr_overrides=pcr_overrides_dict,
            qc_overrides=qc_overrides_dict,
            fixed_prime=req.fixed_prime,         
            mismatch_pos=req.mismatch_pos,       
            mismatch_intensity=req.mismatch_intensity 
        )
        
        ## 🔥 Step 3: 라우터 최적화 (불필요한 재조립 제거)
        ## 이미 스크립트 단에서 완벽한 포맷으로 래핑해 주므로 그대로 반환하되, 시간과 UI 변수만 추가 기록합니다.
        if raw_result.get("status") == "success":
            
            raw_result.setdefault("metadata", {})
            raw_result["metadata"]["timestamp"] = datetime.now().isoformat(timespec="seconds")
            
            raw_result["inputs"] = {
                "pcr_params": pcr_overrides_dict,
                "qc_criteria": qc_overrides_dict,
                "aspcr_settings": {
                    "fixed_prime": req.fixed_prime,
                    "mismatch_pos": req.mismatch_pos,
                    "mismatch_intensity": req.mismatch_intensity
                }
            }
            
            passed = raw_result.get("summary", {}).get("passed_count", 0)
            total = raw_result.get("summary", {}).get("total_count", 0)
            print(f"✅ AS-PCR Design Success. Passed {passed} / {total}")
        else:
            error_msg = raw_result.get("summary", {}).get("error_msg", raw_result.get("reason", "Unknown error"))
            print(f"❌ AS-PCR Design Failed: {error_msg}")
            
        return raw_result

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))