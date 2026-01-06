# app/main.py

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from app.routers import export, qc, pages, design_qpcr, design_methyl, design_aspcr, design_manual

app = FastAPI(
    title="GCX - Primer Design API",
    description="qPCRdesigner 기반 primer 설계 웹 API",
    version="0.1.0",
)

# static /templates 설정
app.mount("/static", StaticFiles(directory="app/static"), name="static")
templates = Jinja2Templates(directory="app/templates", auto_reload=True)

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# -----------------------
#   라우터 등록
# -----------------------
app.include_router(design_qpcr.router)
app.include_router(design_methyl.router)
app.include_router(design_aspcr.router)
app.include_router(design_manual.router)

app.include_router(export.router)
app.include_router(qc.router)   # 👈 QC 라우터
app.include_router(pages.router)   # 👈 QC 라우터


# -----------------------
#   기본 페이지 (GET /)
# -----------------------
@app.get("/", response_class=HTMLResponse)
async def root(request: Request):
    """
    처음 접속할 때는 Primer Design 페이지를 기본으로 보여줌.
    - mode: 'single' 또는 'multi' (템플릿에서 current_design_mode 로 사용)
    """
    return templates.TemplateResponse(
        "design.html",   # 🔹 이제 index.html 대신 design.html 사용
        {
            "request": request,
            "mode": "single",          # 기본 디자인 모드
            "primer_type": "default",  # 필요 시 템플릿에서 사용
            "reference": "hg19",
            "probe": "no",
            "single_result": None,
            "multi_results": None,
            "error": None,
        },
    )
