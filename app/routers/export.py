# app/routers/export.py

from __future__ import annotations

import json
from io import BytesIO
from datetime import datetime
from typing import Any, Dict, List, Optional

import pandas as pd
from fastapi import APIRouter, Form, HTTPException
from fastapi.responses import StreamingResponse
from openpyxl.styles import Font

# internal service (필요하면 사용)
from app.service.amplicone_normlizer import (
    normalize_record,
    df_to_normalized_records,
)


router = APIRouter(
    prefix="/export",
    tags=["export"],
)


@router.post("/excel")
async def export_excel(
    kind: str = Form(...),
    data_json: str = Form(...),
    meta_json: Optional[str] = Form(None),
    styled: Optional[bool] = Form(False),
):
    """
    /export/excel

    - data_json : Amplicon 리스트 (list[dict]) → 엑셀 'amplicons' 시트
    - meta_json : (옵션) 분석 시점 / 파라미터 정보 → 'meta' 시트
    - styled    : (옵션) True 이면 forward/reverse/probe sequence 컬럼에 색상 표시

    kind:
      - 파일 이름 prefix 및 sheet 이름에 사용
    """
    try:
        # -----------------------
        # 1) data_json → DataFrame
        # -----------------------
        rows: List[Dict[str, Any]] = json.loads(data_json)
        if not isinstance(rows, list):
            raise ValueError("data_json must be a list of dicts")

        df = pd.DataFrame(rows)

        # -----------------------
        # 2) meta_json → dict (선택)
        # -----------------------
        meta = None
        if meta_json:
            meta = json.loads(meta_json)

        # -----------------------
        # 3) 엑셀 작성
        # -----------------------
        buffer = BytesIO()
        with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
            # 3-1) 결과 시트
            sheet_name = kind[:31] if kind else "amplicons"
            df.to_excel(writer, index=False, sheet_name=sheet_name)

            # 3-2) 메타 시트 (옵션)
            if meta:
                flat_items: List[Dict[str, Any]] = []

                def flatten(prefix: str, obj: Any):
                    if isinstance(obj, dict):
                        for k, v in obj.items():
                            new_prefix = f"{prefix}.{k}" if prefix else k
                            flatten(new_prefix, v)
                    else:
                        flat_items.append({"key": prefix, "value": obj})

                flatten("", meta)
                meta_df = pd.DataFrame(flat_items)
                meta_df.to_excel(writer, index=False, sheet_name="meta")

            # 3-3) 스타일링 (옵션)
            if styled:
                wb = writer.book
                ws = writer.sheets[sheet_name]

                # 컬럼 인덱스 찾기 (1-based)
                col_map = {name: idx + 1 for idx, name in enumerate(df.columns)}

                red_font = Font(color="FF0000")
                blue_font = Font(color="0000FF")

                # 1행은 헤더이므로 2부터 데이터 row
                for row_idx in range(2, len(df) + 2):
                    # forward / reverse 셀을 빨간색
                    if "forward_sequence" in col_map:
                        cell = ws.cell(row=row_idx, column=col_map["forward_sequence"])
                        cell.font = red_font
                    if "reverse_sequence" in col_map:
                        cell = ws.cell(row=row_idx, column=col_map["reverse_sequence"])
                        cell.font = red_font

                    # probe 셀을 파란색
                    if "probe_sequence" in col_map:
                        cell = ws.cell(row=row_idx, column=col_map["probe_sequence"])
                        cell.font = blue_font

        buffer.seek(0)
        filename = f"{kind}_result_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"

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