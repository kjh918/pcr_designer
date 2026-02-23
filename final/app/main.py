import subprocess
import json
import os
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# -----------------------------------------------------------------------
# 1. API 엔드포인트 (이게 최상단에 있어야 함)
# -----------------------------------------------------------------------
class QpcrRequest(BaseModel):
    chrom: str
    start: int
    end: int
    ref: str
    alt: str
    fasta: str = "/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa"
    top_k: int = 5

@app.post("/api/design/qpcr")
async def run_design(req: QpcrRequest):
    command = [
        "python3.11", "scripts/design_qpcr.py",
        "--chrom", req.chrom,
        "--start", str(req.start),
        "--end", str(req.end),
        "--ref", req.ref,
        "--alt", req.alt,
        "--fasta", req.fasta,
        "--top_k", str(req.top_k)
    ]
    try:
        env = os.environ.copy()
        env["PYTHONPATH"] = "."
        result = subprocess.run(command, capture_output=True, text=True, env=env, check=True)
        return json.loads(result.stdout)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# -----------------------------------------------------------------------
# 2. 루트 경로 접속 시 HTML 반환 (정적 파일 마운트보다 위에 있어야 함)
# -----------------------------------------------------------------------
@app.get("/")
async def serve_home():
    # 경로를 절대 경로로 확인하여 디버깅 용이하게 설정
    base_path = os.path.join(os.getcwd(), "web", "_site")
    main_html = os.path.join(base_path, "main.html")
    index_html = os.path.join(base_path, "index.html")

    if os.path.exists(main_html):
        return FileResponse(main_html)
    elif os.path.exists(index_html):
        return FileResponse(index_html)
    else:
        return {
            "error": "HTML 파일을 찾을 수 없습니다.",
            "checked_paths": [main_html, index_html],
            "current_dir_files": os.listdir(base_path) if os.path.exists(base_path) else "Directory not found"
        }

# -----------------------------------------------------------------------
# 3. 나머지 정적 파일(CSS, JS, 이미지) 마운트
# -----------------------------------------------------------------------
if os.path.exists("web/_site"):
    # html=True 옵션을 빼고 마운트하여 루트 핸들러와 충돌을 방지합니다.
    app.mount("/", StaticFiles(directory="web/_site"), name="site")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8080)