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


router = APIRouter(
	prefix="/export",
	tags=["export"],
)


@router.post("/excel")
async def export_excel(
	kind: str = Form(...),

	# ✅ 고정: sheet1 / sheet2
	total_json: str = Form(...),
	filtered_json: str = Form(...),

	# ✅ 고정: sheet3
	meta_json: Optional[str] = Form(None),

	styled: Optional[bool] = Form(True),
):
	"""
	/export/excel

	고정 시트 구조:
	  - Sheet1 : total_df
	  - Sheet2 : filtered_df
	  - Sheet3 : meta (옵션)

	total_json / filtered_json:
	  - list[dict] (DataFrame records)
	meta_json:
	  - dict (분석 파라미터, threshold, 입력정보 등)
	"""
	try:
		# -----------------------
		# 1) JSON → DataFrame
		# -----------------------
		total_rows: List[Dict[str, Any]] = json.loads(total_json)
		filtered_rows: List[Dict[str, Any]] = json.loads(filtered_json)

		if not isinstance(total_rows, list) or not isinstance(filtered_rows, list):
			raise ValueError("total_json / filtered_json must be list[dict]")

		total_df = pd.DataFrame(total_rows)
		filtered_df = pd.DataFrame(filtered_rows)
		# -----------------------
		# 2) meta sheet (옵션)
		# -----------------------
		meta_df = None
		if meta_json:
			meta_obj = json.loads(meta_json)

			flat_items: List[Dict[str, Any]] = []

			def flatten(prefix: str, obj: Any):
				if isinstance(obj, dict):
					for k, v in obj.items():
						print(k,v)
						new_prefix = f"{prefix}.{k}" if prefix else k
						flatten(new_prefix, v)
				else:
					flat_items.append({"key": prefix, "value": obj})

			flatten("", meta_obj)
			meta_df = pd.DataFrame(flat_items)

		# -----------------------
		# 3) Excel 작성
		# -----------------------
		buffer = BytesIO()
		with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
			total_df.to_excel(writer, index=False, sheet_name="Total_Primer_candidates")
			filtered_df.to_excel(writer, index=False, sheet_name="Filtered_Primer_Candidates")

			if meta_df is not None:
				meta_df.to_excel(writer, index=False, sheet_name="Metadata")
			# -----------------------
			# 4) 스타일링 (옵션)
			# -----------------------
			if styled:
				red_font = Font(color="FF0000")
				blue_font = Font(color="0000FF")

				for sheet_name, df in {
					"Total_Primer_candidates": total_df,
					"Filtered_Primer_Candidates": filtered_df,
				}.items():
					ws = writer.sheets[sheet_name]
					col_map = {name: idx + 1 for idx, name in enumerate(df.columns)}

					for row_idx in range(2, len(df) + 2):
						# forward / reverse → red
						for col in ("forward_sequence", "forward_seq", "reverse_sequence", "reverse_seq"):
							if col in col_map:
								ws.cell(row=row_idx, column=col_map[col]).font = red_font

						# probe → blue
						for col in ("probe_sequence", "probe_seq"):
							if col in col_map:
								ws.cell(row=row_idx, column=col_map[col]).font = blue_font

		buffer.seek(0)
		filename = f"{kind}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"

		return StreamingResponse(
			buffer,
			media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
			headers={
				"Content-Disposition": f'attachment; filename="{filename}"'
			},
		)

	except Exception as e:
		raise HTTPException(status_code=400, detail=f"엑셀 export 실패: {e}")
