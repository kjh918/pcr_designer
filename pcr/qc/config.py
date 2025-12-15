# primer/qc/config.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from config.settings import settings

_qc = settings.qc_params

# Thermo / dimer threshold
HAIRPIN_MIN_DG = _qc.HAIRPIN_MIN_DG
HOMODIMER_MIN_DG = _qc.HOMODIMER_MIN_DG
HETERODIMER_MIN_DG = _qc.HETERODIMER_MIN_DG
PRIMER_MIN_DIFF_TM = _qc.PRIMER_MIN_DIFF_TM
PRIMER_MAX_DIFF_TM = _qc.PRIMER_MAX_DIFF_TM
PROBE_MAX_DIFF_TM = _qc.PROBE_MAX_DIFF_TM
PROBE_MIN_DIFF_TM = _qc.PROBE_MIN_DIFF_TM

# BLAST / amplicon 관련
BLAST_ROOT = _qc.BLAST_ROOT
BLAST_BIN_DIR = _qc.BLAST_BIN_DIR
BLASTN = str(_qc.BLASTN)
BLASTDBCMD = str(_qc.BLASTDBCMD)

BLAST_IDENTITY_THRESHOLD = _qc.BLAST_IDENTITY_THRESHOLD
BLAST_MAX_ALIGNMENTS = _qc.BLAST_MAX_ALIGNMENTS
BLAST_LENGTH_THRESHOLD = _qc.BLAST_LENGTH_THRESHOLD
MIN_AMP_BP = _qc.MIN_AMP_BP
MAX_AMP_BP = _qc.MAX_AMP_BP

BLAST_HIT_MAX = getattr(_qc, "BLAST_HIT_MAX", 1)


@dataclass
class QCThresholds:
    PRIMER_MAX_DIFF_TM: float = PRIMER_MAX_DIFF_TM
    PRIMER_MIN_DIFF_TM: float = PRIMER_MIN_DIFF_TM
    PROBE_MAX_DIFF_TM: float = PROBE_MAX_DIFF_TM
    PROBE_MIN_DIFF_TM: float = PROBE_MIN_DIFF_TM
    hairpin_min_dg: float = HAIRPIN_MIN_DG
    homodimer_min_dg: float = HOMODIMER_MIN_DG
    heterodimer_min_dg: float = HETERODIMER_MIN_DG


@dataclass
class BlastQCConfig:
    identity_threshold: float = BLAST_IDENTITY_THRESHOLD
    length_threshold: int = BLAST_LENGTH_THRESHOLD
    min_amp_bp: int = MIN_AMP_BP
    max_amp_bp: int = MAX_AMP_BP
    max_hits: int = BLAST_HIT_MAX
    max_alignments: int = BLAST_MAX_ALIGNMENTS
