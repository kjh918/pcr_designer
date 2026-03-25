import os
import traceback
from datetime import datetime
from typing import Dict, Any

from fastapi import APIRouter, HTTPException

from pcr.designers.base.schema import QCEvalInput
from scripts.design_qc import evaluate_qc_pipeline 

router = APIRouter(
    prefix="/api/design",
    tags=["qc"]
)

@router.post("/qc")
async def design_qc_api(req: QCEvalInput):
    print(f"\n🚀 [API] Sequence QC Evaluation Request: {req.project_name}")
    
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__)) 
    APP_DIR = os.path.dirname(CURRENT_DIR) 
    ROOT_DIR = os.path.dirname(APP_DIR) 

    SYSTEM_YAML_PATH = os.path.join(ROOT_DIR, "pcr", "config", "system.yaml")
    DESIGNER_YAML_PATH = os.path.join(ROOT_DIR, "pcr", "designers", "base", "config.yaml")

    try:
        req_qc = req.qc_criteria
        qc_core_overrides = {
            "hairpin_min_dg": req_qc.get("hairpin_min_dg", -5.0),
            "homodimer_min_dg": req_qc.get("homodimer_min_dg", -6.0),
            "heterodimer_min_dg": req_qc.get("heterodimer_min_dg", -6.0),
            "min_identity": req_qc.get("min_identity", 90.0),
            "blast_identity_threshold": req_qc.get("min_identity", 90.0),
            "min_hit_length": req_qc.get("min_hit_length", 13),
            "blast_max_alignments": req_qc.get("blast_max_alignments", 50),
            "max_alignments": req_qc.get("blast_max_alignments", 50),
            "min_amp_size": req_qc.get("min_amp_size", 50),
            "min_amp_len": req_qc.get("min_amp_size", 50),
            "max_amp_size": req_qc.get("max_amp_size", 300),
            "max_amp_len": req_qc.get("max_amp_size", 300),
            "use_ispcr_check": req_qc.get("use_ispcr_check", False),
            
            # 🔥 프라이머 오버라이드
            "primer": { 
                "min_diff_tm": req_qc.get("primer", {}).get("min_diff_tm", 0.0),
                "max_diff_tm": req_qc.get("primer", {}).get("max_diff_tm", 3.0) 
            },
            
            # 🔥 프로브 오버라이드 (선택하신 Canvas 코드와 변수명을 정확히 일치시킴)
            "probe": {
                "min_primer_probe_tm_diff": req_qc.get("probe", {}).get("min_primer_probe_tm_diff", 5.0),
                "max_primer_probe_tm_diff": req_qc.get("probe", {}).get("max_primer_probe_tm_diff", 10.0),
                "max_probe_poly_g": req_qc.get("probe", {}).get("max_probe_poly_g", 3),
                "avoid_5_prime_g": req_qc.get("probe", {}).get("avoid_5_prime_g", True)
            }
        }

        # 스크립트 실행 (Base Schema 형태 반환)
        raw_result = evaluate_qc_pipeline(
            project_name=req.project_name,
            fwd_seq=req.sequences.forward,
            rev_seq=req.sequences.reverse,
            probe_seq=req.sequences.probe,
            template_seq=req.sequences.template,
            genome=req.reference_genome,
            qc_overrides=qc_core_overrides,
            base_yaml=DESIGNER_YAML_PATH,
            system_yaml=SYSTEM_YAML_PATH
        )
        
        # 🔥 [JSON 다이어트] BLAST 객체 내부에 중복된 무거운 텍스트/서열 블록 제거
        results_list = raw_result.get("results", [])
        for item in results_list:
            blast_data = item.get("qc_details", {}).get("blast", {})
            if not blast_data: continue
            
        final_output = {
            "status": raw_result.get("status", "error"),
            "metadata": {
                "project_name": req.project_name,
                "reference_genome": req.reference_genome,
                "timestamp": datetime.now().isoformat(timespec="seconds")
            },
            "summary": raw_result.get("summary", {}),
            "inputs": {
                "sequences": req.sequences.model_dump(),
                "qc_criteria": req.qc_criteria
            },
            "results": results_list
        }
        if final_output["status"] == "success":
            print(f"✅ QC Success. Passed {final_output['summary'].get('passed_count', 0)} / {final_output['summary'].get('total_count', 0)}")
        return final_output

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))