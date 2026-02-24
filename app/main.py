import os
import json
import traceback
from typing import Optional, Dict, Any

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel

# 💡 선생님의 내부 라이브러리 임포트
from pcr.seq.fetch import GenomicRegion
from pcr.pipelines.qpcr import run_qpcr
from pcr.config.runtime import (
    get_fasta_handle,
    build_pcr_params_from_web,
    build_qc_params_from_web,
    resolve_pcr_params,
    merge_dict,
)

app = FastAPI()

# CORS 설정
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# ---------------------------------------------------------
# 1. Genome 설정 (경로 매핑)
# ---------------------------------------------------------
GENOME_DB = {
    "hg38": {
        "fasta": "/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa",
        # 필요시 blast 경로 추가
    },
    "hg19": {
        "fasta": "/storage/references_and_index/hg19/fasta/hg19.fa",
    }
}

# ---------------------------------------------------------
# 2. 데이터 모델 (JS에서 보내는 모든 필드 정의)
# ---------------------------------------------------------
class QpcrRequest(BaseModel):
    # --- 기본 좌표 정보 ---
    reference: str = "hg38"
    chrom: str
    start: int
    end: int
    ref: str
    alt: str
    strand: str = "+"
    top_k: int = 5

    # --- Amplicon Size ---
    min_amplicon_length: int = 80
    max_amplicon_length: int = 200

    # --- Primer3 Options (Primer) ---
    primer_min_length: int = 20
    primer_opt_length: int = 25
    primer_max_length: int = 30
    
    primer_min_tm: float = 55.0
    primer_opt_tm: float = 60.0
    primer_max_tm: float = 65.0
    
    primer_min_gc: float = 35.0
    primer_opt_gc: float = 50.0
    primer_max_gc: float = 65.0

    # --- Primer3 Options (Probe) ---
    probe_min_length: int = 20
    probe_opt_length: int = 25
    probe_max_length: int = 30
    
    probe_min_tm: float = 60.0
    probe_opt_tm: float = 65.0 # 보통 Probe Tm은 Primer보다 높게 잡음
    probe_max_tm: float = 70.0
    
    probe_min_gc: float = 35.0
    probe_opt_gc: float = 50.0
    probe_max_gc: float = 65.0
    
    # --- QC / Tm Diff ---
    min_primer_probe_tm_diff: float = 5.0
    max_primer_probe_tm_diff: float = 10.0


# ---------------------------------------------------------
# 3. API 엔드포인트
# ---------------------------------------------------------
@app.post("/api/design/qpcr")
async def design_qpcr_api(req: QpcrRequest):
    print("\n" + "="*50)
    print(f"🚀 [API] qPCR Design 요청 수신: {req.chrom}:{req.start}-{req.end} ({req.strand})")
    print(f"🧬 Genome: {req.reference}")
    print("="*50 + "\n")

    try:
        # 1. Genome 경로 확인
        if req.reference not in GENOME_DB:
            raise HTTPException(status_code=400, detail=f"Unsupported Genome: {req.reference}")
        
        fasta_path = GENOME_DB[req.reference]["fasta"]
        fasta = get_fasta_handle(fasta_path)

        # 2. GenomicRegion 객체 생성
        gr = GenomicRegion(
            chrom=req.chrom, 
            start=req.start, 
            end=req.end, 
            strand=req.strand,
            name=f"{req.reference}_{req.chrom}_{req.start}"
        )

        # -----------------------------------------------------
        # ✅ 3. Config Overrides 구성 (JS 입력값 -> Python Config)
        # -----------------------------------------------------
        # Primer3는 내부적으로 대문자 키(PRIMER_OPT_TM 등)를 사용합니다.
        # 여기서 웹 입력값을 Primer3가 이해하는 구조로 변환합니다.
        
        pcr_overrides = {
            "primer_kwargs": {
                "min_amplicon_length": req.min_amplicon_length,
                "max_amplicon_length": req.max_amplicon_length,
                "n_primers": req.top_k,
                
                # 프라이머 관련 상세 옵션
                "primer3_global_args": {
                    "PRIMER_MIN_SIZE": req.primer_min_length,
                    "PRIMER_OPT_SIZE": req.primer_opt_length,
                    "PRIMER_MAX_SIZE": req.primer_max_length,
                    
                    "PRIMER_MIN_TM": req.primer_min_tm,
                    "PRIMER_OPT_TM": req.primer_opt_tm,
                    "PRIMER_MAX_TM": req.primer_max_tm,
                    
                    "PRIMER_MIN_GC": req.primer_min_gc,
                    "PRIMER_OPT_GC": req.primer_opt_gc,
                    "PRIMER_MAX_GC": req.primer_max_gc,
                }
            },
            
            "probe_kwargs": {
                "n_probes": 1, # TaqMan이므로 기본 1개 이상
                
                # 프로브 관련 상세 옵션 (Primer3에서는 INTERNAL_ 접두어 사용)
                "primer3_global_args": {
                    "PRIMER_INTERNAL_MIN_SIZE": req.probe_min_length,
                    "PRIMER_INTERNAL_OPT_SIZE": req.probe_opt_length,
                    "PRIMER_INTERNAL_MAX_SIZE": req.probe_max_length,
                    
                    "PRIMER_INTERNAL_MIN_TM": req.probe_min_tm,
                    "PRIMER_INTERNAL_OPT_TM": req.probe_opt_tm,
                    "PRIMER_INTERNAL_MAX_TM": req.probe_max_tm,
                    
                    "PRIMER_INTERNAL_MIN_GC": req.probe_min_gc,
                    "PRIMER_INTERNAL_MAX_GC": req.probe_max_gc,
                }
            },
            "bisulfite": {"run": False},
        }

        # 4. QC 파라미터 구성 (Tm 차이 등)
        qc_overrides = {
            "PROBE_MIN_DIFF_TM": req.min_primer_probe_tm_diff,
            "PROBE_MAX_DIFF_TM": req.max_primer_probe_tm_diff,
            "MIN_AMP_BP": req.min_amplicon_length,
            "MAX_AMP_BP": req.max_amplicon_length,
        }

        # 5. 설정 빌드 (Base Config + Overrides)
        # build_pcr_params_from_web 내부에서 기본 yaml을 읽고 pcr_overrides로 덮어씁니다.
        pcr_cfg = build_pcr_params_from_web(pcr_overrides)
        qc_cfg = build_qc_params_from_web(qc_overrides)

        # 6. 최종 파라미터 Resolve (충돌 방지 및 병합)
        resolved = resolve_pcr_params(
            pcr_cfg=pcr_cfg,
            min_amplicon_length=req.min_amplicon_length,
            max_amplicon_length=req.max_amplicon_length,
            n_probes=1,
            n_primers=req.top_k,
            bisulfite=False,
        )
        
        # 7. Primer3 인자 병합
        primer3_global_args = merge_dict(
            base=pcr_cfg.primer_kwargs.primer3_global_args,
            override=None,
        )
        probe_primer3_global_args = merge_dict(
            base=pcr_cfg.probe_kwargs.primer3_global_args,
            override=None,
        )

        # -----------------------------------------------------
        # 🚀 8. 분석 실행
        # -----------------------------------------------------
        result = run_qpcr(
            region=gr,
            fasta=fasta,
            pcr_cfg=pcr_cfg,
            qc_params=qc_cfg,
            min_amplicon_length=resolved.min_amplicon_length,
            max_amplicon_length=resolved.max_amplicon_length,
            n_probes=resolved.n_probes,
            n_primers=resolved.n_primers,
            primer3_global_args=primer3_global_args,
            probe_primer3_global_args=probe_primer3_global_args,
        )

        # 9. 결과 반환 (JSON 리스트)
        return result.total_df.to_dict(orient="records")

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


# ---------------------------------------------------------
# 4. 정적 파일 마운트 (Quarto UI 연결)
# ---------------------------------------------------------
@app.get("/")
async def serve_home():
    # 경로 수정: .web/_site -> web/_site (프로젝트 구조에 맞게)
    base_path = os.path.join(os.getcwd(), "_site") 
    main_html = os.path.join(base_path, "main.html") # 혹은 index.html

    if os.path.exists(main_html):
        return FileResponse(main_html)
    elif os.path.exists(os.path.join(base_path, "index.html")):
        return FileResponse(os.path.join(base_path, "index.html"))
    else:
        return {
            "error": "HTML Not Found", 
            "message": "Quarto 빌드(_site) 폴더를 확인하세요."
        }

if os.path.exists("_site"):
    app.mount("/", StaticFiles(directory="_site"), name="site")

if __name__ == "__main__":
    import uvicorn
    # PYTHONPATH=. python3 main.py 형태로 실행 권장
    uvicorn.run(app, host="0.0.0.0", port=8080)