from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Protocol

import pandas as pd

from pcr.components import Amplicon
from pcr.qc.thermo import evaluate_amplicons
from pcr.config.schema.qc import QCParams


class Designer(Protocol):
    amplicon_list: List[Amplicon]

    def design(self) -> List[Amplicon]: ...
    def reset(self) -> None: ...


@dataclass(frozen=True)
class PipelineResult:
    genomic_id: str
    total_df: pd.DataFrame
    filtered_df: pd.DataFrame


def run_pipeline_from_amplicons(
    *,
    genomic_id: str,
    amplicon_list: List[Amplicon],
    qc_params: QCParams,
    assay: str = "qpcr",
) -> PipelineResult:
    if assay == "qpcr":
        total_rows, filtered_rows = evaluate_amplicons(
            genomic_id,
            amplicon_list,
            qc_params=qc_params,
        )
    else:
        total_rows, filtered_rows = [], []

    total_df = pd.DataFrame(total_rows)
    filtered_df = pd.DataFrame(filtered_rows)

    total_df.index = [genomic_id] * len(total_df)
    filtered_df.index = [genomic_id] * len(filtered_df)

    return PipelineResult(genomic_id, total_df, filtered_df)


def run_pipeline(
        *,
        genomic_id: str,
        qc_params: QCParams,
        assay: str = "qpcr",
        designer: Optional[Designer] = None,
        amplicon_list: Optional[List[Amplicon]] = None,
    ) -> PipelineResult:
    """
    호환용 래퍼:
    - amplicon_list가 주어지면: designer 없이 QC만 수행
    - amplicon_list가 없으면: designer.design()으로 생성 후 QC 수행
    """
    # ✅ 1) amplicon_list가 있으면 그걸 우선 사용
    if amplicon_list is not None and len(amplicon_list) > 0:
        return run_pipeline_from_amplicons(
            genomic_id=genomic_id,
            amplicon_list=amplicon_list,
            qc_params=qc_params,
            assay=assay,
        )

    # ✅ 2) 없으면 designer로부터 생성
    if designer is None:
        raise ValueError("Either 'amplicon_list' must be provided or 'designer' must be provided.")

    designer.design()
    result = run_pipeline_from_amplicons(
        genomic_id=genomic_id,
        amplicon_list=designer.amplicon_list,
        qc_params=qc_params,
        assay=assay,
    )
    designer.reset()
    return result
