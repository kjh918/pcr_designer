# config/schema/qc.py
from __future__ import annotations

from pathlib import Path
from pydantic import BaseModel


class QCParams(BaseModel):
    # -------------------------------------------------
    # Primer / Probe Tm difference QC
    # -------------------------------------------------
    PRIMER_MAX_DIFF_TM: float = 3.0
    PRIMER_MIN_DIFF_TM: float = 0.0
    PROBE_MAX_DIFF_TM: float = 7.0
    PROBE_MIN_DIFF_TM: float = 5.0

    # -------------------------------------------------
    # Thermodynamics QC (primer3)
    # -------------------------------------------------
    HAIRPIN_MIN_DG: float = -5.0
    HOMODIMER_MIN_DG: float = -6.0
    HETERODIMER_MIN_DG: float = -6.0

    # -------------------------------------------------
    # BLAST executable paths
    # -------------------------------------------------
    BLAST_ROOT: Path
    BLAST_BIN_DIR: Path
    BLASTN: Path
    BLASTDBCMD: Path

    # -------------------------------------------------
    # BLAST QC thresholds
    # -------------------------------------------------
    BLAST_IDENTITY_THRESHOLD: float = 80.0
    BLAST_LENGTH_THRESHOLD: int = 10
    BLAST_MAX_ALIGNMENTS: int = 200

    MIN_AMP_BP: int = 50
    MAX_AMP_BP: int = 300
