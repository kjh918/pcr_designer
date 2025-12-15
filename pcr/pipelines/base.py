from __future__ import annotations
from dataclasses import dataclass
from typing import List, Protocol
import pandas as pd

from pcr.qc import QCThresholds, evaluate_amplicons
from pcr.components import Amplicon

class Designer(Protocol):
    amplicon_list: List[Amplicon]
    def design(self) -> List[Amplicon]: ...
    def reset(self) -> None: ...

@dataclass(frozen=True)
class PipelineResult:
    genomic_id: str
    total_df: pd.DataFrame
    filtered_df: pd.DataFrame

def run_pipeline(*, genomic_id: str, designer: Designer, qc_thresholds: QCThresholds) -> PipelineResult:
    designer.design()
    total_rows, filtered_rows = evaluate_amplicons(genomic_id, designer.amplicon_list, qc_thresholds=qc_thresholds)

    total_df = pd.DataFrame(total_rows)
    filtered_df = pd.DataFrame(filtered_rows)
    total_df.index = [genomic_id] * len(total_df)
    filtered_df.index = [genomic_id] * len(filtered_df)

    designer.reset()
    return PipelineResult(genomic_id, total_df, filtered_df)
