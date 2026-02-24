import os
import json
import yaml
import tempfile
import traceback
from typing import Optional, Dict, Any

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

# ✅ Schema Import
from pcr.config.schema.pcr import PCRParams, PrimerKwargs, ProbeKwargs
# ✅ Script Function Import
from scripts.design_qpcr import design_qpcr_primers 

app = FastAPI()

# ----------------------------------------------------------------
# 1. CORS 설정
# ----------------------------------------------------------------
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# ----------------------------------------------------------------
# 2. Request Data Model
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

    # --- Amplicon (Base YAML 값에 맞춤) ---
    min_amplicon_length: int = 60   # [FIX] 80 -> 60
    max_amplicon_length: int = 150  # [FIX] 200 -> 150

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
    
    # [Diff 값] (YAML 생성 시 제거됨)
    min_primer_probe_tm_diff: float = 5.0
    max_primer_probe_tm_diff: float = 10.0
    
    # ✅ [FIX] Base YAML 값에 맞춤 (Absolute Tm)
    probe_min_tm: float = 65.0  # [FIX] 55.0 -> 65.0
    probe_opt_tm: float = 67.0  # [FIX] 60.0 -> 67.0
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
# 3. API 엔드포인트
# ----------------------------------------------------------------
@app.post("/api/design/qpcr")
async def design_qpcr_api(req: QpcrRequest):
    print(f"\n🚀 [API] qPCR Design Request: {req.chrom}:{req.start} ({req.ref}>{req.alt})")
    
    temp_config_path = None
    
    # ----------------------------------------------------------------
    # [경로 설정] 요청하신 대로 dirname 2번 사용
    # ----------------------------------------------------------------
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    SYSTEM_YAML_PATH = os.path.join(BASE_DIR, "pcr", "config", "system.yaml")
    BASE_QPCR_PATH = os.path.join(BASE_DIR, "pcr", "config", "base_pcr.yaml")
    
    print(f"📂 Base Path: {BASE_DIR}")
    print(f"📂 Config Path: {BASE_QPCR_PATH}")

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
            # Diff 값 (삭제될 예정)
            "min_tm_diff": req.min_primer_probe_tm_diff,
            "max_tm_diff": req.max_primer_probe_tm_diff,
            "opt_tm_diff": (req.min_primer_probe_tm_diff + req.max_primer_probe_tm_diff) / 2,
            # Absolute Tm 값 (이게 남음)
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

        # Step 3: Pydantic Schema 검증 및 객체 생성
        pcr_params = PCRParams(
            primer_kwargs=PrimerKwargs(**primer_input_data),
            probe_kwargs=ProbeKwargs(**probe_input_data)
        )
        
        # ----------------------------------------------------------------
        # [핵심 로직] YAML 구조 맞추기 (Diff 제거)
        # ----------------------------------------------------------------
        pcr_params_dict = pcr_params.model_dump()
        
        if pcr_params_dict.get("probe_kwargs"):
            pk = pcr_params_dict["probe_kwargs"]
            # Diff 관련 키 삭제 -> YAML에 절대값만 남게 됨 (base_pcr.yaml 구조 준수)
            pk.pop("min_tm_diff", None)
            pk.pop("max_tm_diff", None)
            pk.pop("opt_tm_diff", None)
            
            print(f"🔧 [DEBUG] Config Sanitized: Using Absolute Tm ({pk.get('min_tm')}-{pk.get('max_tm')})")

        # Step 4: Config Injection (Overwrite)
        if "assay_overrides" not in full_config: full_config["assay_overrides"] = {}
        if "qpcr" not in full_config["assay_overrides"]: full_config["assay_overrides"]["qpcr"] = {}
        
        # 정제된 dict 주입
        full_config["assay_overrides"]["qpcr"]["pcr_params"] = pcr_params_dict

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

        # Step 5: 임시 설정 파일 저장 (⚠️ 삭제 안함)
        with tempfile.NamedTemporaryFile(mode='w', suffix=".yaml", delete=False) as tmp:
            yaml.dump(full_config, tmp)
            temp_config_path = tmp.name

        print(f"📄 [KEEP] Generated Config: {temp_config_path}")

        # Step 6: 디자인 스크립트 실행 (Padding 300 유지)
        raw_result = design_qpcr_primers(
            chrom=req.chrom,
            start=req.start,
            end=req.end,
            ref=req.ref,
            alt=req.alt,
            strand=req.strand,
            genome=req.reference,
            fasta_path=None,
            padding=300, 
            top_k=req.top_k,
            base_yaml=BASE_QPCR_PATH,
            system_yaml=SYSTEM_YAML_PATH 
        )
        
        # Step 7: 결과 반환
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
                }
                formatted_results.append(flat_record)
            
            raw_result["results"] = formatted_results
            print(f"✅ Design Success: Found {len(formatted_results)} candidates.")

        else:
            print(f"❌ Design Failed: {raw_result.get('reason')}")
            if "log_messages" in raw_result:
                print(raw_result["log_messages"])

        return raw_result

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))
    
    finally:
        # ⚠️ Temp 파일 보존 (주석 처리됨)
        if temp_config_path:
             print(f"🔒 Keeping temp config for debugging: {temp_config_path}")
             # os.remove(temp_config_path)

# ----------------------------------------------------------------
# 4. 정적 파일 마운트 (Quarto _site 연결)
# ----------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STATIC_DIR = os.path.join(BASE_DIR, "web", "_site")

print("\n" + "="*50)
print(f"📂 [DEBUG] Frontend Directory: {STATIC_DIR}")

if os.path.exists(STATIC_DIR) and os.path.isdir(STATIC_DIR):
    print(f"✅ FOUND! Serving files from here.")
    app.mount("/", StaticFiles(directory=STATIC_DIR, html=True), name="site")
else:
    print(f"❌ NOT FOUND! Please run 'quarto render' inside the 'web' folder.")
    @app.get("/")
    def index_fallback():
        return {"error": "Frontend not built", "message": f"Path '{STATIC_DIR}' does not exist."}
    
print("="*50 + "\n")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8080)