import os
import yaml
import tempfile
import traceback
from datetime import datetime
from typing import Optional, Dict, Any, List

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from pcr.config.schema.pcr import PCRParams, PrimerKwargs, ProbeKwargs
# 앞서 작성한 매뉴얼 디자인 함수를 import 합니다.
from scripts.validate_primer import design_manual_qpcr 

router = APIRouter(
    prefix="/api/design/manual",
    tags=["manual"]
)
# routers/design_manual.py

class ManualRequest(BaseModel):
    # --- Basic Info ---
    design_name: str
    raw_sequence: str
    target_start: int
    target_end: int
    reference: str = "None"
    top_k: int = 5

    # --- Amplicon ---
    min_amplicon_length: int = 60
    max_amplicon_length: int = 150

    # --- Primer Options ---
    primer_min_length: int = 20
    primer_opt_length: int = 25
    primer_max_length: int = 30
    primer_min_tm: float = 55.0
    primer_opt_tm: float = 60.0
    primer_max_tm: float = 65.0
    primer_min_gc: float = 35.0
    primer_opt_gc: float = 50.0
    primer_max_gc: float = 65.0

    # --- Probe Options ---
    probe_min_length: int = 20
    probe_opt_length: int = 25
    probe_max_length: int = 30
    
    min_primer_probe_tm_diff: float = 5.0
    max_primer_probe_tm_diff: float = 10.0
    
    probe_min_tm: float = 65.0
    probe_opt_tm: float = 67.0
    probe_max_tm: float = 70.0
    
    probe_min_gc: float = 35.0
    probe_opt_gc: float = 50.0
    probe_max_gc: float = 65.0

    # Constraints
    probe_max_poly_g: int = 3
    probe_max_3_end_gc: int = 2
    probe_avoid_5_prime_g: bool = True

    # --- QC Criteria ---
    qc_hairpin_min_dg: float = -5.0
    qc_homodimer_min_dg: float = -6.0
    qc_heterodimer_min_dg: float = -6.0
    qc_min_identity: float = 90.0
    qc_min_hit_length: int = 13
    qc_blast_max_alignments: int = 50
    qc_min_amp_size: int = 50
    qc_max_amp_size: int = 300
    qc_use_ispcr_check: bool = False
    qc_primer_max_diff_tm: float = 3.0
    qc_probe_avoid_5_prime_g: bool = True
    qc_probe_max_poly_g: int = 3
    qc_probe_max_3_end_gc: int = 2

# ----------------------------------------------------------------
# Manual Request Data Model
# ----------------------------------------------------------------
class ManualQpcrRequest(BaseModel):
    # --- Manual Basic Info ---
    design_name: str = "Manual_Project"
    raw_sequence: str  # 사용자가 입력한 전체 DNA 서열
    target_start: int  # 서열 내 프로브 타겟 시작점 (1-based)
    target_end: int    # 서열 내 프로브 타겟 종료점 (1-based)
    reference: str = "none" # BLAST 검증용 레퍼런스
    top_k: int = 5

    # --- Amplicon & Primer ---
    min_amplicon_length: int = 60
    max_amplicon_length: int = 150
    primer_min_length: int = 20
    primer_opt_length: int = 25
    primer_max_length: int = 30
    primer_min_tm: float = 55.0
    primer_opt_tm: float = 60.0
    primer_max_tm: float = 65.0
    primer_min_gc: float = 35.0
    primer_opt_gc: float = 50.0
    primer_max_gc: float = 65.0

    # --- Probe Options ---
    probe_min_length: int = 20
    probe_opt_length: int = 25
    probe_max_length: int = 30
    probe_min_tm: float = 65.0
    probe_opt_tm: float = 67.0
    probe_max_tm: float = 70.0
    probe_min_gc: float = 35.0
    probe_opt_gc: float = 50.0
    probe_max_gc: float = 65.0

    # --- QC & Constraints ---
    probe_max_poly_g: int = 3
    probe_max_3_end_gc: int = 2
    probe_avoid_5_prime_g: bool = True
    qc_hairpin_min_dg: float = -5.0
    qc_homodimer_min_dg: float = -6.0
    qc_heterodimer_min_dg: float = -6.0
    qc_min_identity: float = 90.0
    qc_min_hit_length: int = 13
    qc_blast_max_alignments: int = 50
    qc_use_ispcr_check: bool = False

@router.post("/qpcr")
async def design_manual_api(req: ManualQpcrRequest):
    print(f"🚀 [API] Manual qPCR Design: {req.design_name} (Seq Length: {len(req.raw_sequence)})")
    
    temp_config_path = None
    
    # 경로 설정
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
    ROOT_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR)) 
    SYSTEM_YAML_PATH = os.path.join(ROOT_DIR, "pcr", "config", "system.yaml")
    BASE_QPCR_PATH = os.path.join(ROOT_DIR, "pcr", "config", "base_pcr.yaml")

    try:
        # Step 1: Config Override 설정 (기존 로직과 동일)
        with open(BASE_QPCR_PATH, 'r') as f:
            full_config = yaml.safe_load(f)

        primer_input = {
            "min_amplicon_length": req.min_amplicon_length,
            "max_amplicon_length": req.max_amplicon_length,
            "min_length": req.primer_min_length,
            "opt_length": req.primer_opt_length,
            "max_length": req.primer_max_length,
            "min_tm": req.primer_min_tm,
            "opt_tm": req.primer_opt_tm,
            "max_tm": req.primer_max_tm,
            "min_gc": req.primer_min_gc,
            "opt_gc": req.primer_opt_gc,
            "max_gc": req.primer_max_gc,
        }

        probe_input = {
            "min_length": req.probe_min_length,
            "opt_length": req.probe_opt_length,
            "max_length": req.probe_max_length,
            "min_tm": req.probe_min_tm,
            "opt_tm": req.probe_opt_tm,
            "max_tm": req.probe_max_tm,
            "min_gc": req.probe_min_gc,
            "opt_gc": req.probe_opt_gc,
            "max_gc": req.probe_max_gc,
            "max_probe_poly_g": req.probe_max_poly_g,
            "max_probe_3_end_gc": req.probe_max_3_end_gc,
            "avoid_5_prime_g": req.probe_avoid_5_prime_g
        }

        # Config 주입
        if "assay_overrides" not in full_config: full_config["assay_overrides"] = {}
        full_config["assay_overrides"]["qpcr"] = {
            "pcr_params": {
                "primer_kwargs": primer_input,
                "probe_kwargs": probe_input
            }
        }
        
        full_config["qc_criteria"] = {
            "hairpin_min_dg": req.qc_hairpin_min_dg,
            "homodimer_min_dg": req.qc_homodimer_min_dg,
            "heterodimer_min_dg": req.qc_heterodimer_min_dg,
            "min_identity": req.qc_min_identity,
            "min_hit_length": req.qc_min_hit_length,
            "blast_max_alignments": req.qc_blast_max_alignments,
            "use_ispcr_check": req.qc_use_ispcr_check
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix=".yaml", delete=False) as tmp:
            yaml.dump(full_config, tmp)
            temp_config_path = tmp.name

        # Step 2: Manual Design 실행
        raw_result = design_manual_qpcr(
            design_name=req.design_name,
            raw_sequence=req.raw_sequence,
            target_start=req.target_start,
            target_end=req.target_end,
            genome='none',
            top_k=req.top_k,
            base_yaml=temp_config_path,
            system_yaml=SYSTEM_YAML_PATH
        )

        # Step 3: 결과 데이터 포맷팅
        if raw_result.get("status") == "success":
            formatted_results = []
            for item in raw_result.get("results", []):
                oligos = item.get("oligos", {})
                fwd, rev, prb = oligos.get("forward", {}), oligos.get("reverse", {}), oligos.get("probe", {})
                
                formatted_results.append({
                    "rank": item.get("rank"),
                    "id": item.get("id"),
                    "forward_primer": fwd.get("sequence", "-"),
                    "reverse_primer": rev.get("sequence", "-"),
                    "probe": prb.get("sequence", "-"),
                    "amplicon_size": item.get("amplicon_info", {}).get("length", 0),
                    "tm_f": fwd.get("tm", "-"),
                    "tm_r": rev.get("tm", "-"),
                    "tm_p": prb.get("tm", "-"),
                    "gc_f": fwd.get("gc", "-"),
                    "gc_r": rev.get("gc", "-"),
                    "gc_p": prb.get("gc", "-"),
                    "penalty": item.get("metrics", {}).get("pair_penalty", 0),
                    "alignment_text_block": item.get("alignment_text_block", "")
                })

            return {
                "status": "success",
                "single_result": {
                    "design_name": req.design_name,
                    "total_count": len(formatted_results)
                },
                "single_total_amplicons": formatted_results,
                "export_meta": {
                    "assay": "manual_qpcr",
                    "timestamp": datetime.now().isoformat(timespec="seconds"),
                    "reference_for_blast": req.reference,
                    "target_range": f"{req.target_start}-{req.target_end}",
                    "pcr_params": full_config["assay_overrides"]["qpcr"]["pcr_params"],
                    "qc_params": full_config["qc_criteria"]
                }
            }
        else:
            return {"status": "fail", "error": raw_result.get("reason"), "log_messages": raw_result.get("log_messages", "")}

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        if temp_config_path and os.path.exists(temp_config_path):
            os.remove(temp_config_path)