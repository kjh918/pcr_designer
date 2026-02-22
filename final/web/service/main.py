import os
import sys

from fastapi import FastAPI, Request, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

# 프로젝트 루트 등록
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from pcr.factory import PCRFactory
from pcr.config.loader import load_pipeline_config

app = FastAPI(title="GENCURIX PCR Web Service")

# CORS (필요시)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# template / static 설정
templates = Jinja2Templates(directory=os.path.join(BASE_DIR, "web/templates"))
app.mount("/static", StaticFiles(directory=os.path.join(BASE_DIR, "web/static")), name="static")


# -----------------------
# UI 페이지
# -----------------------
@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


# -----------------------
# Primer Design API
# -----------------------
@app.post("/api/design/qpcr")
async def design_qpcr(request: Request):

    base_yaml = os.path.join(BASE_DIR, "pcr/config/base_pcr.yaml")
    sys_yaml = os.path.join(BASE_DIR, "pcr/config/system.yaml")

    try:
        web_input = await request.json()

        config = load_pipeline_config(
            base_yaml_path=base_yaml,
            system_yaml_path=sys_yaml,
            assay_type="qPCR"
        )

        factory = PCRFactory(config)

        result = factory.execute_design(
            chrom=web_input["chrom"],
            start=web_input["start"],
            end=web_input["end"]
        )

        return result.model_dump()

    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "web.service.main:app",
        host="0.0.0.0",
        port=8000,
        reload=True
    )