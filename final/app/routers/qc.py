# app/routers/design_qc.py
import os
import yaml
import tempfile
import traceback
from datetime import datetime
from typing import Optional, Dict, Any, List

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from pcr.config.schema.pcr import PCRParams, PrimerKwargs, ProbeKwargs
from scripts.design_qc import evaluate_qc_pipeline 

# 🔹 라우터 객체 생성
router = APIRouter(
    prefix="/api/design",
    tags=["qc"]
)

# ----------------------------------------------------------------
# Request Data Model
# ----------------------------------------------------------------
class QpcrRequest(BaseModel):
    # --- Basic Info ---
    reference: str = "hg38"
    chrom: str
    start: int
    end: int
    ref: str
    alt: str
    strand: str = "+"
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
# API 엔드포인트 (POST /api/design/qc)
# ----------------------------------------------------------------
@router.post("/qc")
async def design_qc_api(req: QpcrRequest):
    print(f"\n🚀 [API] qPCR Design Request: {req.chrom}:{req.start} ({req.ref}>{req.alt})")
    
    temp_config_path = None
    
    # [경로 설정] 
    # 이 파일은 app/routers/ 에 있으므로, 프로젝트 루트는 3단계 상위(../ ../ ../)가 아니라
    # 실행 컨텍스트(run.py 위치)에 따라 다르지만, 보통abspath(__file__) 기준 2단계 위(app -> root)
    # 안전하게 계산:
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__)) # app/routers
    APP_DIR = os.path.dirname(CURRENT_DIR) # app
    ROOT_DIR = os.path.dirname(APP_DIR) # project_root

    SYSTEM_YAML_PATH = os.path.join(ROOT_DIR, "pcr", "config", "system.yaml")
    BASE_QPCR_PATH = os.path.join(ROOT_DIR, "pcr", "config", "base_pcr.yaml")

    try:
        # Step 1: Base YAML 로드
        if not os.path.exists(BASE_QPCR_PATH):
            raise FileNotFoundError(f"Config file not found: {BASE_QPCR_PATH}")
        
        with open(BASE_QPCR_PATH, 'r') as f:
            full_config = yaml.safe_load(f)

        # Step 2: Request -> Dict Mapping
        primer_input_data = {
            "n_candidates": 100,
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

        probe_input_data = {
            "n_candidates": 100,
            "min_length": req.probe_min_length,
            "opt_length": req.probe_opt_length,
            "max_length": req.probe_max_length,
            "min_tm_diff": req.min_primer_probe_tm_diff,
            "max_tm_diff": req.max_primer_probe_tm_diff,
            "opt_tm_diff": (req.min_primer_probe_tm_diff + req.max_primer_probe_tm_diff) / 2,
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

        # Step 3: Pydantic Schema
        pcr_params = PCRParams(
            primer_kwargs=PrimerKwargs(**primer_input_data),
            probe_kwargs=ProbeKwargs(**probe_input_data)
        )
        
        # Step 4: Diff 제거 및 Config 주입
        pcr_params_dict = pcr_params.model_dump()
        if pcr_params_dict.get("probe_kwargs"):
            pk = pcr_params_dict["probe_kwargs"]
            pk.pop("min_tm_diff", None)
            pk.pop("max_tm_diff", None)
            pk.pop("opt_tm_diff", None)
            print(f"🔧 [DEBUG] Config Sanitized: Using Absolute Tm ({pk.get('min_tm')}-{pk.get('max_tm')})")

        if "assay_overrides" not in full_config: full_config["assay_overrides"] = {}
        if "qc" not in full_config["assay_overrides"]: full_config["assay_overrides"]["qc"] = {}
        full_config["assay_overrides"]["qc"]["pcr_params"] = pcr_params_dict

        # QC Injection
        if "qc_criteria" not in full_config: full_config["qc_criteria"] = {}
        qc_conf = full_config["qc_criteria"]
        qc_conf["hairpin_min_dg"] = req.qc_hairpin_min_dg
        qc_conf["homodimer_min_dg"] = req.qc_homodimer_min_dg
        qc_conf["heterodimer_min_dg"] = req.qc_heterodimer_min_dg
        qc_conf["min_identity"] = req.qc_min_identity
        qc_conf["min_hit_length"] = req.qc_min_hit_length
        qc_conf["blast_max_alignments"] = req.qc_blast_max_alignments
        qc_conf["min_amp_size"] = req.qc_min_amp_size
        qc_conf["max_amp_size"] = req.qc_max_amp_size
        qc_conf["use_ispcr_check"] = req.qc_use_ispcr_check
        if "primer" not in qc_conf: qc_conf["primer"] = {}
        qc_conf["primer"]["max_diff_tm"] = req.qc_primer_max_diff_tm
        if "probe" not in qc_conf: qc_conf["probe"] = {}
        qc_conf["probe"]["avoid_5_prime_g"] = req.qc_probe_avoid_5_prime_g
        qc_conf["probe"]["max_probe_poly_g"] = req.qc_probe_max_poly_g
        qc_conf["probe"]["max_probe_3_end_gc"] = req.qc_probe_max_3_end_gc

        # Step 5: 임시 파일 저장
        with tempfile.NamedTemporaryFile(mode='w', suffix=".yaml", delete=False) as tmp:
            yaml.dump(full_config, tmp)
            temp_config_path = tmp.name

        print(f"📄 [KEEP] Generated Config: {temp_config_path}")

        # Step 6: 디자인 스크립트 실행
        raw_result = design_qc_primers(
            chrom=req.chrom,
            start=req.start,
            end=req.end,
            ref=req.ref,
            alt=req.alt,
            strand=req.strand,
            genome=req.reference,
            fasta_path=None,
            padding=150, 
            top_k=req.top_k,
            base_yaml=temp_config_path,
            system_yaml=SYSTEM_YAML_PATH 
        )
        
        # Step 7: 결과 데이터 포맷팅
        if raw_result.get("status") == "success":
            formatted_results = []
            for item in raw_result.get("results", []):
                oligos = item.get("oligos", {})
                fwd = oligos.get("forward", {})
                rev = oligos.get("reverse", {})
                prb = oligos.get("probe") or {}
                amp_info = item.get("amplicon_info", {})
                metrics = item.get("metrics", {})

                flat_record = {
                    "rank": item.get("rank"),
                    "id": item.get("id"),
                    "forward_primer": fwd.get("sequence", "-"),
                    "reverse_primer": rev.get("sequence", "-"),
                    "probe": prb.get("sequence", "-"),
                    "amplicon_size": amp_info.get("length", 0),
                    "genomic_pos": amp_info.get("genomic_pos", "-"),
                    "tm_f": fwd.get("tm", "-"),
                    "tm_r": rev.get("tm", "-"),
                    "tm_p": prb.get("tm", "-"),
                    "gc_f": fwd.get("gc", "-"),
                    "gc_r": rev.get("gc", "-"),
                    "gc_p": prb.get("gc", "-"),
                    "penalty": metrics.get("pair_penalty", 0),
                    "alignment_text_block": item.get("alignment_text_block", "")
                }
                formatted_results.append(flat_record)

            export_meta = {
                "assay": "qc",
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "reference": req.reference,
                "region": {
                    "chrom": req.chrom,
                    "start": req.start,
                    "end": req.end,
                    "name": f"{req.ref}>{req.alt}",
                    "strand": req.strand
                },
                "total_count": len(formatted_results),    
                "filtered_count": len(formatted_results), 
                
                # --- 모든 입력 파라미터 통합 주입 ---
                "pcr_params": {
                    "amplicon": {
                        "min_length": req.min_amplicon_length,
                        "max_length": req.max_amplicon_length
                    },
                    "primer": {
                        "min_length": req.primer_min_length,
                        "opt_length": req.primer_opt_length,
                        "max_length": req.primer_max_length,
                        "min_tm": req.primer_min_tm,
                        "opt_tm": req.primer_opt_tm,
                        "max_tm": req.primer_max_tm,
                        "min_gc": req.primer_min_gc,
                        "opt_gc": req.primer_opt_gc,
                        "max_gc": req.primer_max_gc
                    },
                    "probe": {
                        "min_length": req.probe_min_length,
                        "opt_length": req.probe_opt_length,
                        "max_length": req.probe_max_length,
                        "min_gc": req.probe_min_gc,
                        "opt_gc": req.probe_opt_gc,
                        "max_gc": req.probe_max_gc,
                        "min_tm_diff_from_primer": req.min_primer_probe_tm_diff,
                        "max_tm_diff_from_primer": req.max_primer_probe_tm_diff,
                    }
                },
                "qc_params": {
                    "thermodynamics": {
                        "hairpin_min_dg": req.qc_hairpin_min_dg,
                        "homodimer_min_dg": req.qc_homodimer_min_dg,
                        "heterodimer_min_dg": req.qc_heterodimer_min_dg
                    },
                    "specificity": {
                        "min_identity": req.qc_min_identity,
                        "min_hit_length": req.qc_min_hit_length,
                        "blast_max_alignments": req.qc_blast_max_alignments,
                        "use_ispcr_check": req.qc_use_ispcr_check
                    },
                    "advanced": {
                        "primer_max_diff_tm": req.qc_primer_max_diff_tm,
                        "probe_avoid_5_prime_g": req.qc_probe_avoid_5_prime_g,
                        "probe_max_poly_g": req.qc_probe_max_poly_g,
                        "probe_max_3_end_gc": req.qc_probe_max_3_end_gc
                    }
                }
            }

            response_data = {
                "status": "success",
                "single_result": {
                    "region": export_meta["region"],
                    "total_count": export_meta["total_count"],
                    "filtered_count": export_meta["filtered_count"],
                },
                "single_total_amplicons": formatted_results, 
                "single_filtered_amplicons": formatted_results,
                "export_meta": export_meta
            }
            
            print(f"✅ Design Success: Found {len(formatted_results)} candidates.")
            return response_data

        else:
            print(f"❌ Design Failed: {raw_result.get('reason')}")
            return {
                "status": "fail",
                "error": raw_result.get("reason"),
                "log_messages": raw_result.get("log_messages", ""),
                "single_result": None,
                "single_total_amplicons": [],
                "single_filtered_amplicons": []
            }

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))
    
    finally:
        if temp_config_path:
             print(f"🔒 Keeping temp config: {temp_config_path}")