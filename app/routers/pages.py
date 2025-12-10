from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

templates = Jinja2Templates("app/templates")

router = APIRouter(tags=["pages"])

@router.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse(
        "design.html",
        {
            "request": request,
            "mode": "single",
            "single_result": None,
            "multi_results": None,
            "error": None,
        },
    )

@router.get("/design", response_class=HTMLResponse)
async def design_page(request: Request):
    return templates.TemplateResponse(
        "design.html",
        {
            "request": request,
            "mode": "single",
            "single_result": None,
            "multi_results": None,
            "error": None,
        },
    )

@router.get("/qc", response_class=HTMLResponse)
async def qc_page(request: Request):
    return templates.TemplateResponse(
        "qc.html",
        {
            "request": request,
            "qc_result": None,
            "error": None,
        },
    )