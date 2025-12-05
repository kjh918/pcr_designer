# app/routers/export.py

from __future__ import annotations

import json
from io import BytesIO

import pandas as pd
from fastapi import APIRouter, Form, HTTPException
from fastapi.responses import StreamingResponse

router = APIRouter(
    prefix="/export",
    tags=["export"],
)


# templates = Jinja2Templates(directory="app/templates")

@router.post("/excel")
async def export_amplicons_to_excel(
    kind: str = Form(...),      # "single_filtered", "single_total", "multi_1" 등
    data_json: str = Form(...), # 템플릿에서 넘겨주는 JSON 문자열
):
    """
    템플릿에서 넘겨준 amplicon 리스트(JSON)를 엑셀로 변환하여 다운로드.
    """
    try:
        rows = json.loads(data_json)
        if not isinstance(rows, list):
            raise ValueError("data_json must be a list of dicts")

        df = pd.DataFrame(rows)

        buffer = BytesIO()
        with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
            sheet_name = kind[:31] if kind else "Sheet1"
            df.to_excel(writer, index=False, sheet_name=sheet_name)

        buffer.seek(0)
        filename = f"{kind}_amplicons.xlsx"

        return StreamingResponse(
            buffer,
            media_type=(
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            ),
            headers={
                "Content-Disposition": f'attachment; filename="{filename}"'
            },
        )
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"엑셀 export 실패: {e}")
