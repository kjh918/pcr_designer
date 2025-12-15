# primer/qc/qc_eval.py
from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Tuple

from pcr.qc.config import QCThresholds, BlastQCConfig
from pcr.qc.thermo import evaluate_amplicons
from pcr.qc.blast import blast_qc_for_primer_pair


def run_thermo_qc(
    genomic_id,
    amplicons: Iterable[Any],
    *,
    thresholds: QCThresholds = QCThresholds(),
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    return evaluate_amplicons(genomic_id, amplicons, thresholds)


def run_blast_qc(
    f_name: str,
    f_seq: str,
    r_name: str,
    r_seq: str,
    db: str,
    *,
    probe_name: Optional[str] = None,
    probe_seq: Optional[str] = None,
    config: BlastQCConfig = BlastQCConfig(),
) -> Dict[str, Any]:
    return blast_qc_for_primer_pair(
        f_name, f_seq, r_name, r_seq, db,
        probe_name=probe_name, probe_seq=probe_seq,
        config=config,
    )
