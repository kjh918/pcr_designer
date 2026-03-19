import os
import traceback
from datetime import datetime
from typing import Dict, Any

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

# 🔥 [핵심 수정] 수정 중이신 파일명(run_qc)을 정확히 지정하여 Import 합니다!
from scripts.design_qc import evaluate_qc_pipeline 

router = APIRouter(
    prefix="/api/design",
    tags=["qc"]
)

class SequencesModel(BaseModel):
    forward: str
    reverse: str
    probe: str = ""
    template: str = ""

class QcThermodynamicsModel(BaseModel):
    hairpin_min_dg: float = -5.0
    homodimer_min_dg: float = -6.0
    heterodimer_min_dg: float = -6.0

class QcBlastModel(BaseModel):
    min_identity: float = 90.0
    min_hit_length: int = 13
    max_alignments: int = 50

class QcAmpliconModel(BaseModel):
    min_size: int = 50
    max_size: int = 300
    use_ispcr: bool = False

class QcOligoModel(BaseModel):
    max_tm_diff: float = 3.0

class QcCriteriaModel(BaseModel):
    thermodynamics: QcThermodynamicsModel
    blast: QcBlastModel
    amplicon: QcAmpliconModel
    oligo: QcOligoModel

class QCRequest(BaseModel):
    project_name: str = "QC_Project"
    sequences: SequencesModel
    reference_genome: str = "none"
    qc_criteria: QcCriteriaModel


@router.post("/qc")
async def design_qc_api(req: QCRequest):
    print(f"\n🚀 [API] Sequence QC Evaluation Request: {req.project_name}")
    
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__)) 
    APP_DIR = os.path.dirname(CURRENT_DIR) 
    ROOT_DIR = os.path.dirname(APP_DIR) 

    SYSTEM_YAML_PATH = os.path.join(ROOT_DIR, "pcr", "config", "system.yaml")
    BASE_QPCR_PATH = os.path.join(ROOT_DIR, "pcr", "config", "base_pcr.yaml")

    try:
        qc_flat_overrides = {
            "hairpin_min_dg": req.qc_criteria.thermodynamics.hairpin_min_dg,
            "homodimer_min_dg": req.qc_criteria.thermodynamics.homodimer_min_dg,
            "heterodimer_min_dg": req.qc_criteria.thermodynamics.heterodimer_min_dg,
            
            "min_identity": req.qc_criteria.blast.min_identity,
            "min_identity_threshold": req.qc_criteria.blast.min_identity,
            "min_hit_length": req.qc_criteria.blast.min_hit_length,
            "blast_max_alignments": req.qc_criteria.blast.max_alignments,
            
            "min_amp_size": req.qc_criteria.amplicon.min_size,
            "max_amp_size": req.qc_criteria.amplicon.max_size,
            "use_ispcr_check": req.qc_criteria.amplicon.use_ispcr,
            
            "oligo": {"max_tm_diff": req.qc_criteria.oligo.max_tm_diff}
        }

        # 함수 호출
        raw_result = evaluate_qc_pipeline(
            project_name=req.project_name,
            fwd_seq=req.sequences.forward,
            rev_seq=req.sequences.reverse,
            probe_seq=req.sequences.probe,
            template_seq=req.sequences.template,
            genome=req.reference_genome,
            qc_overrides=qc_flat_overrides,
            base_yaml=BASE_QPCR_PATH,
            system_yaml=SYSTEM_YAML_PATH
        )
        
        # QC 탈락건도 무조건 화면에 띄우기 위해 필터링 스킵
        if "single_total_amplicons" in raw_result:
            raw_result["single_filtered_amplicons"] = raw_result["single_total_amplicons"]

        if raw_result.get("status") == "success":
            export_meta = {
                "assay": "qc",
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "project_name": req.project_name,
                "reference": req.reference_genome,
                "total_count": 1,    
                "filtered_count": len(raw_result.get("single_filtered_amplicons", [])), 
                "sequences": req.sequences.model_dump(),
                "qc_params": qc_flat_overrides
            }

            raw_result["export_meta"] = export_meta
            
            print(f"✅ QC Evaluation Success. Extracted {len(raw_result['single_filtered_amplicons'])} blocks.")
            return raw_result

        else:
            print(f"❌ QC Evaluation Failed: {raw_result.get('reason')}")
            return raw_result

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))