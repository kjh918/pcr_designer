# pcr/qc/runner.py
from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Tuple

from pcr.config.schema.qc import QCParams  # 타입/값 객체만
from pcr.qc.thermo import evaluate_amplicons
from pcr.qc.blast import blast_qc_for_primer_pair


def run_thermo_qc(
	genomic_id: str,
	amplicons: Iterable[Any],
	*,
	qc_params: QCParams,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
	return evaluate_amplicons(genomic_id, amplicons, qc_params=qc_params)


def run_blast_qc(
	f_name: str,
	f_seq: str,
	r_name: str,
	r_seq: str,
	db: str,
	fasta: str,
	*,
	qc_params: QCParams,
	probe_name: Optional[str] = None,
	probe_seq: Optional[str] = None,
) -> Dict[str, Any]:
	return blast_qc_for_primer_pair(
		f_name=f_name,
		f_seq=f_seq,
		r_name=r_name,
		r_seq=r_seq,
		db=db,
		fasta=fasta,
		qc_params=qc_params,
		probe_name=probe_name,
		probe_seq=probe_seq,
	)
