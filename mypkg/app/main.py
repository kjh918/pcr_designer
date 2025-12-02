# app/main.py

from fastapi import FastAPI, Request, Form, UploadFile, File, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles

from app.schemas import (
    SingleRegionDesignRequest,
    SingleRegionDesignResponse,
    MultiRegionDesignResult,
    RegionInput,
)
from primer.core import design_qpcr_for_region  # ★ qPCR wrapper
# from primer.core import design_qpcr_for_region  # ★ qPCR wrapper

import pandas as pd
from io import BytesIO


app = FastAPI(
    title="GCX - Primer Design API",
    description="qPCRdesigner 기반 primer 설계 웹 API",
    version="0.1.0",
)

# static / templates 설정
app.mount("/static", StaticFiles(directory="app/static"), name="static")
templates = Jinja2Templates(directory="app/templates")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# -----------------------
#   기본 페이지 (GET /)
# -----------------------
@app.get("/", response_class=HTMLResponse)
async def root(request: Request):
    """처음 접속할 때는 결과 없이 빈 폼만 보여줌"""
    return templates.TemplateResponse(
        "index.html",
        {
            "request": request,
            "single_result": None,
            "single_primers": None,   # ★ 처음에는 비움
            "multi_results": None,
            "error": None,
        },
    )

# -------------------------------
#   /design/form (GET) → / 리다이렉트
# -------------------------------
@app.get("/design/form", response_class=HTMLResponse)
async def design_form_get():
    """주소창에 /design/form 을 직접 치면 / 로 돌려보냄"""
    return RedirectResponse(url="/")


@app.post("/design/form", response_class=HTMLResponse)
async def design_from_form(
    request: Request,
    mode: str = Form("single"),           # single / multi
    primer_type: str = Form("default"),   # 지금은 안 써도 일단 유지
    reference: str = Form("hg19"),        # hg19 / hg38
    probe: str = Form("no"),              # yes / no  ← 기본값 No (웹 기본 상태와 맞춤)

    # --- Single region inputs ---
    chrom: str = Form(""),
    start: int | None = Form(None),
    end: int | None = Form(None),
    name: str = Form(""),

    # --- Primer 설정 (None이면 config 기본값 사용) ---
    min_amplicon_length: int | None = Form(None),
    max_amplicon_length: int | None = Form(None),
    n_primers: int | None = Form(None),

    primer_opt_length: int | None = Form(None),
    primer_min_length: int | None = Form(None),
    primer_max_length: int | None = Form(None),
    primer_opt_gc: float | None = Form(None),
    primer_min_gc: float | None = Form(None),
    primer_max_gc: float | None = Form(None),

    # --- Probe 설정 (None이면 config 기본값 / 또는 미사용) ---
    n_probes: int | None = Form(None),
    min_primer_probe_tm_diff: float | None = Form(None),
    max_primer_probe_tm_diff: float | None = Form(None),
    probe_opt_length: int | None = Form(None),
    probe_min_length: int | None = Form(None),
    probe_max_length: int | None = Form(None),
    probe_opt_tm: float | None = Form(None),
    probe_min_tm: float | None = Form(None),
    probe_max_tm: float | None = Form(None),
    probe_opt_gc: float | None = Form(None),
    probe_min_gc: float | None = Form(None),
    probe_max_gc: float | None = Form(None),

    # --- Multi 모드용 파일 ---
    file: UploadFile | None = File(None),
):
    single_result = None
    single_primers = None   # ★ 테이블용 후보 리스트
    multi_results = None
    error = None

    try:
        # probe 토글에 따라 실제 사용할 n_probes 결정
        # - probe == "no" 면 n_probes를 강제로 0으로 덮어써서 probe 비활성화
        if probe == "yes":
            effective_n_probes = n_probes  # None이면 core에서 config 값 사용
        else:
            effective_n_probes = 0

        # =======================
        #   1) 단일 영역 모드
        # =======================
        if mode == "single":
            if not chrom or start is None or end is None:
                raise HTTPException(
                    status_code=400,
                    detail="qPCR 설계에는 chrom / start / end 가 모두 필요합니다.",
                )

            # RegionInput은 1-based 좌표를 쓰도록 유지
            region = RegionInput(
                chrom=chrom,
                start=start,
                end=end,
                name=name or None,
                sequence="",       # qPCR에서는 사용하지 않음
            )

            # ★ best_pair + 후보 리스트 동시 반환
            if probe == 'yes':

                primer_pair, primer_rows = design_qpcr_for_region(
                    region=region,
                    reference_name=reference,

                    # ---- Primer high-level ----
                    min_amplicon_length=min_amplicon_length,
                    max_amplicon_length=max_amplicon_length,
                    n_primers=n_primers,

                    # ---- Probe high-level ----
                    n_probes=effective_n_probes,

                    # ---- Primer 세부 옵션 (None이면 config 값 사용) ----
                    primer_opt_length=primer_opt_length,
                    primer_min_length=primer_min_length,
                    primer_max_length=primer_max_length,
                    primer_opt_gc=primer_opt_gc,
                    primer_min_gc=primer_min_gc,
                    primer_max_gc=primer_max_gc,

                    # ---- Probe 세부 옵션 (None이면 config 값 사용) ----
                    min_primer_probe_tm_diff=min_primer_probe_tm_diff,
                    max_primer_probe_tm_diff=max_primer_probe_tm_diff,
                    probe_opt_length=probe_opt_length,
                    probe_min_length=probe_min_length,
                    probe_max_length=probe_max_length,
                    probe_opt_tm=probe_opt_tm,
                    probe_min_tm=probe_min_tm,
                    probe_max_tm=probe_max_tm,
                    probe_opt_gc=probe_opt_gc,
                    probe_min_gc=probe_min_gc,
                    probe_max_gc=probe_max_gc,
                )

                single_result = SingleRegionDesignResponse(
                    region=region,
                    primer_pair=primer_pair,
                )
                single_primers = primer_rows   # ★ 템플릿으로 넘길 리스트
            else:
                
                primer_pair, primer_rows = design_qpcr_for_region(
                    region=region,
                    reference_name=reference,

                    # ---- Primer high-level ----
                    min_amplicon_length=min_amplicon_length,
                    max_amplicon_length=max_amplicon_length,
                    n_primers=n_primers,

                    # ---- Probe high-level ----
                    n_probes=effective_n_probes,

                    # ---- Primer 세부 옵션 (None이면 config 값 사용) ----
                    primer_opt_length=primer_opt_length,
                    primer_min_length=primer_min_length,
                    primer_max_length=primer_max_length,
                    primer_opt_gc=primer_opt_gc,
                    primer_min_gc=primer_min_gc,
                    primer_max_gc=primer_max_gc,

                    # ---- Probe 세부 옵션 (None이면 config 값 사용) ----
                    probe_opt_length=probe_opt_length,
                    probe_min_length=probe_min_length,
                    probe_max_length=probe_max_length,
                    probe_opt_tm=probe_opt_tm,
                    probe_min_tm=probe_min_tm,
                    probe_max_tm=probe_max_tm,
                    probe_opt_gc=probe_opt_gc,
                    probe_min_gc=probe_min_gc,
                    probe_max_gc=probe_max_gc,
                )

        # =======================
        #   2) 다중(Excel) 모드
        # =======================
        elif mode == "multi":
            if file is None or not file.filename:
                raise HTTPException(
                    status_code=400,
                    detail="Multiple 모드에서는 Excel 파일이 필요합니다.",
                )

            contents = await file.read()
            df = pd.read_excel(BytesIO(contents))

            required_cols = ["chrom", "start", "end"]
            for col in required_cols:
                if col not in df.columns:
                    raise HTTPException(
                        status_code=400,
                        detail=f"필수 컬럼이 없습니다: {col}",
                    )

            multi_results: list[MultiRegionDesignResult] = []

            for idx, row in df.iterrows():
                chrom_val = str(row["chrom"])
                start_val = int(row["start"])
                end_val = int(row["end"])
                name_val = (
                    str(row["name"])
                    if "name" in df.columns and not pd.isna(row["name"])
                    else None
                )

                region = RegionInput(
                    chrom=chrom_val,
                    start=start_val,   # 엑셀도 1-based 로 들어온다고 가정
                    end=end_val,
                    name=name_val,
                    sequence="",   # qPCR에서는 사용하지 않음
                )

                primer_pair, _ = design_qpcr_for_region(
                    region=region,
                    reference_name=reference,

                    min_amplicon_length=min_amplicon_length,
                    max_amplicon_length=max_amplicon_length,
                    n_primers=n_primers,
                    n_probes=effective_n_probes,

                    primer_opt_length=primer_opt_length,
                    primer_min_length=primer_min_length,
                    primer_max_length=primer_max_length,
                    primer_opt_gc=primer_opt_gc,
                    primer_min_gc=primer_min_gc,
                    primer_max_gc=primer_max_gc,

                    min_primer_probe_tm_diff=min_primer_probe_tm_diff,
                    max_primer_probe_tm_diff=max_primer_probe_tm_diff,
                    probe_opt_length=probe_opt_length,
                    probe_min_length=probe_min_length,
                    probe_max_length=probe_max_length,
                    probe_opt_tm=probe_opt_tm,
                    probe_min_tm=probe_min_tm,
                    probe_max_tm=probe_max_tm,
                    probe_opt_gc=probe_opt_gc,
                    probe_min_gc=probe_min_gc,
                    probe_max_gc=probe_max_gc,
                )

                multi_results.append(
                    MultiRegionDesignResult(
                        region=region,
                        primer_pair=primer_pair,
                    )
                )

        else:
            raise HTTPException(
                status_code=400,
                detail=f"알 수 없는 mode: {mode}",
            )

    except ValueError as e:
        error = f"Reference/설계 에러: {e}"
    except HTTPException as e:
        error = e.detail
    except Exception as e:
        error = str(e)

    return templates.TemplateResponse(
        "index.html",
        {
            "request": request,
            "single_result": single_result,
            "single_primers": single_primers,
            "multi_results": multi_results,
            "error": error,
        },
    )
