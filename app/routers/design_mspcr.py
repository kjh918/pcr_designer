import os
import traceback
from datetime import datetime
from fastapi import APIRouter, HTTPException

# 🔥 전용 API 스키마 임포트
from pcr.designers.ms_pcr.schema import MSPCRDesignInput
# 코어 디자인 스크립트 임포트
from scripts.design_mspcr import design_mspcr_primers 

router = APIRouter(
    prefix="/api/design",
    tags=["mspcr"]
)

@router.post("/mspcr")
async def design_mspcr_api(req: MSPCRDesignInput):
    """
    MS-PCR 디자인 엔드포인트.
    """
    print(f"\n🚀 [API] MS-PCR Design Request: {req.design_name}")
    
    # [시스템 경로 설정]
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__)) 
    APP_DIR = os.path.dirname(CURRENT_DIR) 
    ROOT_DIR = os.path.dirname(APP_DIR) 

    SYSTEM_YAML_PATH = os.path.join(ROOT_DIR, "pcr", "config", "system.yaml")
    BASE_PCR_PATH = os.path.join(ROOT_DIR, "pcr", "config", "base_pcr.yaml")

    try:
        # 🔥 Step 1: UI에서 받은 모든 파라미터(Tm, GC, 길이 등)를 딕셔너리로 추출
        pcr_overrides_dict = req.to_core_pcr_params()
        qc_overrides_dict = req.to_core_qc_overrides()
        print(qc_overrides_dict)

        # 🔥 Step 2: 디자인 스크립트 실행
        raw_result = design_mspcr_primers(
            design_name=req.design_name,
            raw_sequence_with_brackets=req.sequence,
            # [수정] 'none' 하드코딩을 제거하고 UI에서 선택한 값을 전달합니다.
            genome=req.reference_genome, 
            top_k=req.top_k,
            window_size_3prime=req.qc_criteria.oligo.window_size_3prime,
            min_cpg_count=req.qc_criteria.oligo.min_cpg_count,
            base_yaml=BASE_PCR_PATH,
            system_yaml=SYSTEM_YAML_PATH,
            # 🔥 UI 파라미터 덮어쓰기(Overwrite)를 위한 인자 전달
            pcr_overrides=pcr_overrides_dict,
            qc_overrides=qc_overrides_dict
        )
        # Step 3: 결과 반환
        if raw_result.get("status") == "success":
            # 결과 메타데이터 보강
            raw_result.setdefault("metadata", {})
            raw_result["metadata"]["project_name"] = req.design_name
            raw_result["metadata"]["reference_genome"] = req.reference_genome
            
            # 내보내기용 메타데이터 추가
            raw_result["export_meta"] = {
                "assay": "mspcr",
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "project_name": req.design_name,
                "reference": req.reference_genome,
                "pcr_params": pcr_overrides_dict,
                "qc_params": qc_overrides_dict
            }
            return raw_result
        else:
            return {"status": "fail", "reason": raw_result.get("reason", "Design failed")}

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))