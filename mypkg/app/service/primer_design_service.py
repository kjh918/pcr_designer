# app/services/primer_design_service.py

from __future__ import annotations

from io import BytesIO
from typing import Any, Dict, List

import pandas as pd

from app.schemas import RegionInput
from primer.core import design_qpcr_for_region
from app.services.amplicon_normalizer import df_to_normalized_records


class PrimerDesignService:
    """
    qPCR primer/probe 설계 서비스.
    - core.design_qpcr_for_region 호출
    - DataFrame → list[dict] 변환 및 정규화
    """

    def design_single_region(
        self,
        region: RegionInput,
        reference_name: str,
        design_kwargs: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        단일 리전 설계.
        """
        total_df, filtered_df = design_qpcr_for_region(
            region=region,
            reference_name=reference_name,
            **design_kwargs,
        )

        total_records = df_to_normalized_records(total_df)
        filtered_records = df_to_normalized_records(filtered_df)

        return {
            "region": region,
            "total_amplicons": total_records,
            "filtered_amplicons": filtered_records,
            "total_count": len(total_records),
            "filtered_count": len(filtered_records),
        }

    def design_multi_from_excel(
        self,
        excel_bytes: bytes,
        reference_name: str,
        design_kwargs: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """
        엑셀 파일 (bytes)을 받아 multi-region 설계 수행.
        엑셀은 chrom / start / end / (name) 컬럼을 가진다고 가정.
        """
        df = pd.read_excel(BytesIO(excel_bytes))

        required_cols = ["chrom", "start", "end"]
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"필수 컬럼이 없습니다: {col}")

        results: List[Dict[str, Any]] = []

        for _, row in df.iterrows():
            chrom_val = str(row["chrom"])
            start_val = int(row["start"])
            end_val = int(row["end"])
            name_val = (
                str(row["name"])
                if "name" in df.columns and not pd.isna(row["name"])
                else None
            )

            region = RegionInput(
                chrom=chrom_val,
                start=start_val,
                end=end_val,
                name=name_val,
                sequence="",
            )

            total_df, filtered_df = design_qpcr_for_region(
                region=region,
                reference_name=reference_name,
                **design_kwargs,
            )

            results.append(
                {
                    "region": region,
                    "total_amplicons": df_to_normalized_records(total_df),
                    "filtered_amplicons": df_to_normalized_records(filtered_df),
                    "total_count": len(total_df),
                    "filtered_count": len(filtered_df),
                }
            )

        return results
