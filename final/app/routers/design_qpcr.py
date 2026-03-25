import os
import traceback
from datetime import datetime
from fastapi import APIRouter, HTTPException

# 🔥 최신 Pydantic 기반 API Request 스키마 임포트
from pcr.designers.qpcr.schema import QPCRDesignInput
from scripts.design_qpcr import design_qpcr_pipeline 

router = APIRouter(
    prefix="/api/design",
    tags=["qpcr"]
)

@router.post("/qpcr")
async def design_qpcr_api(req: QPCRDesignInput):
    print(f"\n🚀 [API] qPCR Design Request: {req.design_name}")
    
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__)) 
    APP_DIR = os.path.dirname(CURRENT_DIR) 
    ROOT_DIR = os.path.dirname(APP_DIR) 

    # 시스템 및 qPCR 전용 config 경로 지정
    SYSTEM_YAML_PATH = os.path.join(ROOT_DIR, "pcr", "config", "system.yaml")
    BASE_QPCR_PATH = os.path.join(ROOT_DIR, "pcr", "designers", "qpcr", "config.yaml")

    try:
        # 🔥 지저분한 매핑, Temp YAML 생성 로직을 모두 지우고 Schema 객체 하나만 던집니다.
        raw_result = design_qpcr_pipeline(
            req=req,
            designer_yaml=BASE_QPCR_PATH,
            system_yaml=SYSTEM_YAML_PATH
        )
        
        if raw_result.get("status") == "success":
            # 라우터에서 공통 응답 규격의 메타데이터 보강
            raw_result.setdefault("metadata", {})
            raw_result["metadata"]["project_name"] = req.design_name
            raw_result["metadata"]["reference_genome"] = req.reference_genome
            raw_result["metadata"]["timestamp"] = datetime.now().isoformat(timespec="seconds")
            
            # summary 블록에서 통계 추출 (프론트엔드 출력을 위함)
            passed = raw_result.get("summary", {}).get("passed_count", 0)
            total = raw_result.get("summary", {}).get("total_count", 0)
            
            if passed > 0:
                print(f"✅ [Execution Success] qPCR Design Passed: {passed} / {total} candidates.")
            else:
                print(f"⚠️ [Execution Success] qPCR Design Completed, but ALL FAILED: {passed} / {total} candidates.")
        else:
            print(f"❌ [Execution Failed] qPCR Pipeline Error: {raw_result.get('reason')}")

        return raw_result

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))