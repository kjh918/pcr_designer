# app/services/amplicon_normalizer.py

from __future__ import annotations
from typing import Any, Dict, List
import pandas as pd

STANDARD_COLUMNS = [
    "reference_template_sequence",
    "amplicon_sequence",
    "amplicon_length",

    "forward_sequence",
    "reverse_sequence",
    "probe_sequence",         # probe 없으면 "" 로 채움

    "forward_tm",
    "reverse_tm",
    "probe_tm",

    "forward_gc_percent",
    "reverse_gc_percent",
    "probe_gc_percent",

    "forward_hairpin_tm",
    "reverse_hairpin_tm",
    "probe_hairpin_tm",

    "forward_homodimer_dg",
    "reverse_homodimer_dg",
    "probe_homodimer_dg",

    "heterodimer_dg",
    "heterodimer_tm",
]


def normalize_record(record: Dict[str, Any]) -> Dict[str, Any]:
    """
    Amplicon.to_dict() 결과 하나(dict)를 받아
    STANDARD_COLUMNS 기준으로 key가 없으면 ""로 채워 넣는다.
    """
    normalized: Dict[str, Any] = {}
    for col in STANDARD_COLUMNS:
        normalized[col] = record.get(col, "")
    return normalized


def df_to_normalized_records(df: pd.DataFrame) -> List[Dict[str, Any]]:
    """
    DataFrame → list[dict] → 각 dict를 normalize_record로 정규화.
    """
    records = df.to_dict(orient="records")
    return [normalize_record(r) for r in records]
