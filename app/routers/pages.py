from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

templates = Jinja2Templates("app/templates")

router = APIRouter(tags=["pages"])


def _normalize_assay(assay: str) -> str:
    if assay in ("qpcr", "methyl", "as-pcr"):
        return assay
    return "qpcr"


def _assay_to_primer_type(assay: str) -> str:
    if assay == "methyl":
        return "methyl"
    if assay == "as-pcr":
        return "as"
    return "default"


@router.get("/", response_class=HTMLResponse, name="home_page")
async def home(request: Request):
    # 홈은 qpcr 기본으로 design 페이지 보여주는 형태 유지
    assay = "qpcr"
    return templates.TemplateResponse(
        "design.html",
        {
            "request": request,
            "mode": "single",
            "assay": assay,
            "primer_type": _assay_to_primer_type(assay),
            "reference": "hg19",
            "probe": "no",
            "single_result": None,
            "multi_results": None,
            "error": None,
        },
    )


@router.get("/design", response_class=HTMLResponse, name="design_page")
async def design_page(request: Request, assay: str = "qpcr"):
    assay = _normalize_assay(assay)

    return templates.TemplateResponse(
        "design.html",
        {
            "request": request,
            "mode": "single",
            "assay": assay,  # ✅ 이걸로 qpcr/methyl/as-pcr 분기 & active 표시 가능
            "primer_type": _assay_to_primer_type(assay),
            "reference": "hg19",
            "probe": "no",
            "single_result": None,
            "multi_results": None,
            "error": None,
        },
    )


@router.get("/qc", response_class=HTMLResponse, name="qc_page")
async def qc_page(request: Request):
    return templates.TemplateResponse(
        "qc.html",
        {
            "request": request,
            "qc_result": None,
            "error": None,
        },
    )
