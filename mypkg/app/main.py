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


# --------------------------------------------
#   /design/form (POST)  : HTML 폼 처리 (qPCR 전용)
# --------------------------------------------
@app.post("/design/form", response_class=HTMLResponse)
async def design_from_form(
    request: Request,
    mode: str = Form("single"),           # single / multi
    primer_type: str = Form("default"),   # 지금은 안 써도 일단 유지
    reference: str = Form("hg19"),        # hg19 / hg38 (sequence 모드는 qPCR에서 지원 X)
    probe: str = Form("yes"),             # yes / no (qPCRdesigner 내부에서 probe 사용)
    chrom: str = Form(""),
    start: int | None = Form(None),
    end: int | None = Form(None),
    name: str = Form(""),

    min_amplicon_length: int = Form(80),
    max_amplicon_length: int = Form(120),

    file: UploadFile | None = File(None),
):
    single_result = None
    single_primers = None   # ★ 테이블용 후보 리스트
    multi_results = None
    error = None

    try:
        # qPCRdesigner 는 reference FASTA 기반이라 sequence 모드는 지원 X
        if reference == "sequence":
            raise HTTPException(
                status_code=400,
                detail="qPCR 모드에서는 reference=sequence 는 지원하지 않습니다. hg19/hg38 등의 reference genome 을 선택하세요.",
            )

        # =======================
        #   1) 단일 영역 모드
        # =======================
        if mode == "single":
            if not chrom or start is None or end is None:
                raise HTTPException(
                    status_code=400,
                    detail="qPCR 설계에는 chrom / start / end 가 모두 필요합니다.",
                )
            

            region = RegionInput(
                chrom=chrom,
                start=start-1,
                end=end,
                name=name or None,
                sequence="",   # qPCR에서는 사용하지 않음
            )
            print(region)
            # ★ best_pair + 후보 리스트 동시 반환
            primer_pair, primer_rows = design_qpcr_for_region(
                region,
                reference_name=reference,
                min_amplicon_length=min_amplicon_length,
                max_amplicon_length=max_amplicon_length,
                # 필요하면 n_probes, n_primers, bisulfite, n_best 등도 추가
            )

            single_result = SingleRegionDesignResponse(
                region=region,
                primer_pair=primer_pair,
            )
            single_primers = primer_rows   # ★ 템플릿으로 넘길 리스트

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
                    start=start_val,
                    end=end_val,
                    name=name_val,
                    sequence="",   # qPCR에서는 사용하지 않음
                )

                # 다중 모드에서는 일단 최상위 1쌍만 사용
                primer_pair, _ = design_qpcr_for_region(
                    region,
                    reference_name=reference,
                    min_amplicon_length=min_amplicon_length,
                    max_amplicon_length=max_amplicon_length,
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
            "single_primers": single_primers,   # ★ 여기서 넘김
            "multi_results": multi_results,
            "error": error,
        },
    )
