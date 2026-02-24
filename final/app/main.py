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

# ✅ Schema Import (수정하지 않음)
from pcr.config.schema.pcr import PCRParams, PrimerKwargs, ProbeKwargs
# ✅ Script Function Import
from scripts.design_qpcr import design_qpcr_primers 

app = FastAPI()

# 1. CORS 설정
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# 2. Request Data Model
# Frontend(JS)에서 보내는 필드명과 일치해야 합니다.
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
    min_amplicon_length: int = 80
    max_amplicon_length: int = 200

    # --- Primer Options (Frontend prefixes: primer_) ---
    primer_min_length: int = 20
    primer_opt_length: int = 25
    primer_max_length: int = 30
    
    primer_min_tm: float = 55.0
    primer_opt_tm: float = 60.0
    primer_max_tm: float = 65.0
    
    primer_min_gc: float = 35.0
    primer_opt_gc: float = 50.0
    primer_max_gc: float = 65.0

    # --- Probe Options (Frontend prefixes: probe_) ---
    probe_min_length: int = 20
    probe_opt_length: int = 25
    probe_max_length: int = 30
    
    # Probe Tm은 Diff(차이)로 관리되므로 Diff 값을 받습니다.
    min_primer_probe_tm_diff: float = 5.0
    max_primer_probe_tm_diff: float = 10.0
    
    probe_min_gc: float = 35.0
    probe_opt_gc: float = 50.0
    probe_max_gc: float = 65.0


@app.post("/api/design/qpcr")
async def design_qpcr_api(req: QpcrRequest):
    print(f"\n🚀 [API] qPCR Design Request: {req.chrom}:{req.start} ({req.ref}>{req.alt})")
    
    temp_config_path = None
    try:
        # ----------------------------------------------------------------
        # 1. Base YAML 로드
        # ----------------------------------------------------------------

        base_yaml_path = f"{os.getcwd()}/pcr/config/base_pcr.yaml"
        if not os.path.exists(base_yaml_path):
            raise FileNotFoundError(f"Config file not found: {base_yaml_path}")

        with open(base_yaml_path, 'r') as f:
            full_config = yaml.safe_load(f)

        # ----------------------------------------------------------------
        # 2. Request -> Schema Mapping (Main.py에서 처리)
        # ----------------------------------------------------------------
        
        # (1) PrimerKwargs용 데이터 구성
        # Request의 'primer_min_tm' -> Schema의 'min_tm' 으로 매핑
        primer_input_data = {
            "n_candidates": 100, # 내부 기본값
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

        # (2) ProbeKwargs용 데이터 구성
        # Request의 'probe_min_length' -> Schema의 'min_length' 등으로 매핑
        probe_input_data = {
            "n_candidates": 100,
            "min_length": req.probe_min_length,
            "opt_length": req.probe_opt_length,
            "max_length": req.probe_max_length,
            
            # Tm Diff 매핑
            "min_tm_diff": req.min_primer_probe_tm_diff,
            "max_tm_diff": req.max_primer_probe_tm_diff,
            # opt_tm_diff는 min/max 사이의 값으로 자동설정하거나, 입력을 안 받으면 기본값 사용
            "opt_tm_diff": (req.min_primer_probe_tm_diff + req.max_primer_probe_tm_diff) / 2,

            "min_gc": req.probe_min_gc,
            "opt_gc": req.probe_opt_gc,
            "max_gc": req.probe_max_gc,
        }

        # ----------------------------------------------------------------
        # 3. Pydantic Schema 검증 및 객체 생성
        # ----------------------------------------------------------------
        # 여기서 값이 유효하지 않으면(예: 문자열 들어옴) 에러가 발생하여 안전함
        pcr_params = PCRParams(
            primer_kwargs=PrimerKwargs(**primer_input_data),
            probe_kwargs=ProbeKwargs(**probe_input_data)
        )

        # ----------------------------------------------------------------
        # 4. YAML 구조 업데이트 (Injection)
        # ----------------------------------------------------------------
        # assay_overrides -> qpcr -> pcr_params 섹션 덮어쓰기
        if "assay_overrides" not in full_config:
            full_config["assay_overrides"] = {}
        if "qpcr" not in full_config["assay_overrides"]:
            full_config["assay_overrides"]["qpcr"] = {}
            
        # model_dump()를 사용해 깔끔한 dict로 변환 후 주입
        full_config["assay_overrides"]["qpcr"]["pcr_params"] = pcr_params.model_dump()

        # ----------------------------------------------------------------
        # 5. 임시 설정 파일 저장
        # ----------------------------------------------------------------
        with tempfile.NamedTemporaryFile(mode='w', suffix=".yaml", delete=False) as tmp:
            yaml.dump(full_config, tmp)
            temp_config_path = tmp.name

        print(f"📄 Generated Temp Config: {temp_config_path}")

        # ----------------------------------------------------------------
        # 6. 디자인 스크립트 실행
        # ----------------------------------------------------------------
        raw_result = design_qpcr_primers(
            chrom=req.chrom,
            start=req.start,
            end=req.end,
            ref=req.ref,
            alt=req.alt,
            strand=req.strand,
            genome=req.reference,  # 예: hg38
            fasta_path=None,       # Config(system.yaml)에서 자동 조회
            padding=100,
            top_k=req.top_k,
            base_yaml=temp_config_path,  # ✅ 생성한 임시 Config 경로 전달
            system_yaml="pcr/config/system.yaml"
        )
        
        # ----------------------------------------------------------------
        # 7. 결과 데이터 변환 (Flattening for Frontend)
        # ----------------------------------------------------------------
        if raw_result.get("status") == "success":
            formatted_results = []
            
            for item in raw_result.get("results", []):
                # 중첩 구조 분해
                oligos = item.get("oligos", {})
                fwd = oligos.get("forward", {})
                rev = oligos.get("reverse", {})
                prb = oligos.get("probe") or {}
                
                amp_info = item.get("amplicon_info", {})
                metrics = item.get("metrics", {})

                # Frontend Table용 Flat Data 생성
                flat_record = {
                    "rank": item.get("rank"),
                    "id": item.get("id"),
                    
                    # Sequences
                    "forward_primer": fwd.get("sequence", "-"),
                    "reverse_primer": rev.get("sequence", "-"),
                    "probe": prb.get("sequence", "-"),
                    
                    # Amplicon Info
                    "amplicon_size": amp_info.get("length", 0),
                    "genomic_pos": amp_info.get("genomic_pos", "-"),
                    
                    # Tm
                    "tm_f": fwd.get("tm", "-"),
                    "tm_r": rev.get("tm", "-"),
                    "tm_p": prb.get("tm", "-"),
                    
                    # GC (Optional)
                    "gc_f": fwd.get("gc", "-"),
                    "gc_r": rev.get("gc", "-"),
                    "gc_p": prb.get("gc", "-"),
                    
                    # Penalty
                    "penalty": metrics.get("pair_penalty", 0),
                }
                formatted_results.append(flat_record)
            
            # 원본 결과 교체
            raw_result["results"] = formatted_results
            print(f"✅ Design Success: Found {len(formatted_results)} candidates.")

        else:
            print(f"❌ Design Failed: {raw_result.get('reason')}")

        return raw_result

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))
    
    finally:
        # 8. 임시 파일 삭제
        if temp_config_path and os.path.exists(temp_config_path):
            os.remove(temp_config_path)
            print("🗑️ Temp config cleaned up.")

# 정적 파일 서빙 (Quarto 빌드 결과물)
if os.path.exists("_site"):
    app.mount("/", StaticFiles(directory="_site", html=True), name="site")
elif os.path.exists("../_site"):
    app.mount("/", StaticFiles(directory="../_site", html=True), name="site")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8080)