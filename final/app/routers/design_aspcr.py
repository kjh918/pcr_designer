# app/routers/design_aspcr.py
import os
import yaml
import tempfile
import traceback
from datetime import datetime
from typing import Optional, Dict, Any, List

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from pcr.config.schema.pcr import PCRParams, PrimerKwargs
from scripts.design_aspcr import design_aspcr_primers 

# 🔹 라우터 객체 생성
router = APIRouter(
    prefix="/api/design",
    tags=["aspcr"]
)

# ----------------------------------------------------------------
# Request Data Model (AS-PCR 맞춤형으로 수정)
# ----------------------------------------------------------------
class AspcrRequest(BaseModel):
    # --- Basic Info ---
    reference: str = "hg38"
    chrom: str
    start: int
    end: int
    ref: str
    alt: str
    strand: str = "+"
    top_k: int = 5

    # --- AS-PCR Specific Options (🔥 추가됨) ---
    fixed_prime: str = "forward"
    mismatch_pos: int = 3
    mismatch_intensity: str = "strong"

    # --- Amplicon ---
    min_amplicon_length: int = 60
    max_amplicon_length: int = 150

    # --- Primer Options ---
    primer_min_length: int = 15
    primer_opt_length: int = 22
    primer_max_length: int = 30
    primer_min_tm: float = 52.0
    primer_opt_tm: float = 58.0
    primer_max_tm: float = 65.0
    primer_min_gc: float = 35.0
    primer_opt_gc: float = 50.0
    primer_max_gc: float = 65.0

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
    qc_primer_max_diff_tm: float = 5.0 # AS-PCR은 미스매치로 인해 Tm 차이가 크게 날 수 있으므로 넉넉히 줍니다.


# ----------------------------------------------------------------
# API 엔드포인트 (POST /api/design/aspcr)
# ----------------------------------------------------------------
@router.post("/aspcr")
async def design_aspcr_api(req: AspcrRequest):
    print(f"\n🚀 [API] ASPCR Design Request: {req.chrom}:{req.start} ({req.ref}>{req.alt})")
    
    temp_config_path = None
    
    # [경로 설정] 
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

        # Step 3: Pydantic Schema
        pcr_params = PCRParams(
            primer_kwargs=PrimerKwargs(**primer_input_data)
        )
        
        # Step 4: Diff 제거 및 Config 주입 (프로브 관련 코드 싹 삭제됨)
        pcr_params_dict = pcr_params.model_dump()

        if "assay_overrides" not in full_config: full_config["assay_overrides"] = {}
        if "aspcr" not in full_config["assay_overrides"]: full_config["assay_overrides"]["aspcr"] = {}
        full_config["assay_overrides"]["aspcr"]["pcr_params"] = pcr_params_dict

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

        # Step 5: 임시 파일 저장
        with tempfile.NamedTemporaryFile(mode='w', suffix=".yaml", delete=False) as tmp:
            yaml.dump(full_config, tmp)
            temp_config_path = tmp.name

        print(f"📄 [KEEP] Generated Config: {temp_config_path}")

        # Step 6: 디자인 스크립트 실행 (AS-PCR 파라미터 전달)
        raw_result = design_aspcr_primers(
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
            system_yaml=SYSTEM_YAML_PATH,
            fixed_prime=req.fixed_prime,         # 🔥 AS-PCR 변수 전달
            mismatch_pos=req.mismatch_pos,       # 🔥 AS-PCR 변수 전달
            mismatch_intensity=req.mismatch_intensity # 🔥 AS-PCR 변수 전달
        )
        
        # Step 7: 결과 데이터 포맷팅
        if raw_result.get("status") == "success":
            
            # 🔥 [핵심 보정] 우리가 디자인 스크립트에서 정성껏 만든 "세트 단위 데이터 구조"를 
            # 억지로 평탄화(Flatten)하지 않고 그대로 프론트엔드로 전달합니다!
            formatted_sets = raw_result.get("results", [])

            export_meta = {
                "assay": "aspcr",
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "reference": req.reference,
                "region": {
                    "chrom": req.chrom,
                    "start": req.start,
                    "end": req.end,
                    "name": f"{req.ref}>{req.alt}",
                    "strand": req.strand
                },
                "total_count": len(formatted_sets),    
                "filtered_count": len(formatted_sets), 
                
                # --- 모든 입력 파라미터 통합 주입 (프로브 제외됨) ---
                "pcr_params": {
                    "aspcr_settings": {
                        "fixed_prime": req.fixed_prime,
                        "mismatch_pos": req.mismatch_pos,
                        "mismatch_intensity": req.mismatch_intensity
                    },
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
                        "primer_max_diff_tm": req.qc_primer_max_diff_tm
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
                # 프론트엔드 호환성을 위해 키 이름은 유지하되 내용은 Set 구조를 담아 보냅니다.
                "single_total_amplicons": formatted_sets, 
                "single_filtered_amplicons": formatted_sets,
                "export_meta": export_meta
            }
            
            print(f"✅ Design Success: Found {len(formatted_sets)} candidate sets.")
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