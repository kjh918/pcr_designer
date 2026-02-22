import os
import sys
import uvicorn
from fastapi import FastAPI, Request, HTTPException
from fastapi.responses import HTMLResponse # 화면 출력을 위해 추가

# 프로젝트 경로 설정
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE_DIR)

from pcr.factory import PCRFactory
from pcr.config.loader import load_pipeline_config

app = FastAPI(title="GENCURIX PCR API")

# 1. [핵심 수정] 브라우저 접속 시 "Not Found" 대신 "서버 작동 중" 화면을 보여줍니다.
@app.get("/", response_class=HTMLResponse)
async def root():
    return """
    <html>
        <head><title>Gencurix API</title></head>
        <body style="font-family: sans-serif; text-align: center; padding-top: 50px;">
            <h1 style="color: #2b384c;">🚀 Gencurix PCR Backend is Running!</h1>
            <p>API Endpoint: <code>/api/design/qpcr</code> (POST)</p>
            <p>API Documentation: <a href="/docs">Swagger UI (/docs)</a></p>
            <div style="margin-top: 20px; padding: 10px; background: #f8f9fa; display: inline-block; border-radius: 5px;">
                <b>Note:</b> 실제 대시보드는 <code>quarto preview index.qmd</code>로 실행하세요.
            </div>
        </body>
    </html>
    """

# 2. 기존 디자인 로직 (POST 방식)
@app.post("/api/design/qpcr")
async def design_qpcr(request: Request):
    base_yaml = os.path.join(BASE_DIR, "pcr/config/base_pcr.yaml")
    sys_yaml = os.path.join(BASE_DIR, "pcr/config/system.yaml")

    try:
        config = load_pipeline_config(base_yaml, sys_yaml, "qPCR")
        web_input = await request.json()
        
        factory = PCRFactory(config)
        result = factory.execute_design(
            chrom=web_input.get("chrom"),
            start=web_input.get("start"),
            end=web_input.get("end")
        )
        return result.model_dump()
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    # 포트 8000번으로 실행
    uvicorn.run(app, host="0.0.0.0", port=8000)