from pathlib import Path
from pydantic import BaseModel

class QCParams(BaseModel):
    PRIMER_MAX_DIFF_TM: float = 3
    PRIMER_MIN_DIFF_TM: float = 0
    PROBE_MAX_DIFF_TM: float = 7
    PROBE_MIN_DIFF_TM: float = 5

    HAIRPIN_MIN_DG: float = -5.0
    HOMODIMER_MIN_DG: float = -6.0
    HETERODIMER_MIN_DG: float = -6.0

    BLAST_ROOT: Path
    BLAST_BIN_DIR: Path
    BLASTN: Path
    BLASTDBCMD: Path

    BLAST_IDENTITY_THRESHOLD: float = 80.0
    BLAST_MAX_ALIGNMENTS: int = 200
    BLAST_LENGTH_THRESHOLD: int = 10
    MIN_AMP_BP: int = 50
    MAX_AMP_BP: int = 300
