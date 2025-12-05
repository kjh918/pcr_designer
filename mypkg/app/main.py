# app/main.py

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates


from app.routers import design, export


app = FastAPI(
    title="GCX - Primer Design API",
    description="qPCRdesigner 기반 primer 설계 웹 API",
    version="0.1.0",
)

# static /templates 설정
app.mount("/static", StaticFiles(directory="app/static"), name="static")
templates = Jinja2Templates(directory="app/templates")

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# -----------------------
#   /design 관련 라우터
# -----------------------

app.include_router(design.router)
app.include_router(export.router)

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
            "single_primers": None,   # 처음에는 비움
            "multi_results": None,
            "error": None,
        },
    )

