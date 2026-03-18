# app/main.py
from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse
import os

# 라우터 임포트 (파일 분리됨)
from app.routers import design_qpcr
from app.routers import design_manual
from app.routers import design_aspcr
from app.routers import design_mspcr
from app.routers import qc

app = FastAPI(
    title="GCX - Primer Design API",
    description="qPCRdesigner 기반 primer 설계 웹 API",
    version="0.1.0",
)

# -----------------------
#   CORS 설정 (모든 IP 허용)
# -----------------------
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
# design_qpcr.py에 있는 router를 /api/design/qpcr 경로로 연결하는 것이 아니라
# 라우터 내부에서 경로를 정의했으므로 include만 하면 됩니다.
app.include_router(design_qpcr.router)
app.include_router(design_manual.router)
app.include_router(design_aspcr.router)
app.include_router(design_mspcr.router)
app.include_router(qc.router)


# -----------------------
#   정적 파일 (Frontend)
# -----------------------
# 프로젝트 루트 경로 계산
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STATIC_DIR = os.path.join(BASE_DIR, "web", "_site")

if os.path.exists(STATIC_DIR):
    app.mount("/", StaticFiles(directory=STATIC_DIR, html=True), name="static")
else:
    print(f"⚠️ Warning: Frontend build not found at {STATIC_DIR}")

# -----------------------
#   Health Check
# -----------------------
@app.get("/health")
async def health_check():
    return {"status": "ok", "service": "GCX Primer Design API"}