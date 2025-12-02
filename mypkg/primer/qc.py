# primer_qc/thermo.py

import primer3
from Bio.SeqUtils import gc_fraction

from .schema import PrimerThermoResult
# primer_qc/runner.py

import os
from typing import List, Optional

from .thermo import compute_thermo, compute_heterodimer
from .blast import run_blast_for_primers, find_nearby_amplicons
from .qc_rules import compute_qc_flags
from .schema import PrimerThermoResult  # type hint 용도로만 사용


class PrimerThermoBlastQC:
    """
    - 입력 TSV를 읽어서
    - Thermo + BLAST + QC 룰을 적용하고
    - 결과 TSV를 쓰는 오케스트레이터.
    """

    def __init__(
        self,
        input_path: str,
        output_dir: str,
        blast_db: str,
        blastn_path: str = "blastn",
        identity_threshold: float = 85.0,
        length_threshold: int = 12,
        max_alignments: int = 200,
        hairpin_dg_cutoff: float = -5.0,
        homodimer_dg_cutoff: float = -6.0,
        heterodimer_dg_cutoff: float = -6.0,
        tm_min: float = 55.0,
        tm_max: float = 70.0,
        tm_diff_max: float = 3.0,
        blast_hit_max: int = 1,
        min_amp_bp: int = 50,
        max_amp_bp: int = 300,
    ):
        self.input_path = input_path
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.output_file = os.path.join(self.output_dir, "primer_thermo_blast_result.tab")

        self.blast_db = blast_db
        self.blastn_path = blastn_path

        # BLAST 기준
        self.identity_threshold = identity_threshold
        self.length_threshold = length_threshold
        self.max_alignments = max_alignments

        # Thermo 기준
        self.hairpin_dg_cutoff = hairpin_dg_cutoff
        self.homodimer_dg_cutoff = homodimer_dg_cutoff
        self.heterodimer_dg_cutoff = heterodimer_dg_cutoff

        # Tm 기준
        self.tm_min = tm_min
        self.tm_max = tm_max
        self.tm_diff_max = tm_diff_max

        # BLAST hit 개수 기준
        self.blast_hit_max = blast_hit_max

        # 잠재 amplicon 길이 기준
        self.min_amp_bp = min_amp_bp
        self.max_amp_bp = max_amp_bp

    def run(self):
        with open(self.input_path, "r") as handle, open(self.output_file, "w") as out:
            self._write_header(out)

            for line in handle:
                if line.startswith("Forward_Primer"):
                    continue

                row = line.rstrip("\n").split("\t")
                if len(row) < 4:
                    continue

                f_name, f_seq_raw, r_name, r_seq_raw = row[0], row[1], row[2], row[3]
                f_seq = f_seq_raw.strip().upper()
                r_seq = r_seq_raw.strip().upper()

                # '+' 포함된 primer는 스킵
                if "+" in f_seq or "+" in r_seq:
                    continue

                result_line = self._process_primer_pair(f_name, f_seq, r_name, r_seq)
                print(result_line)
                out.write(result_line + "\n")

        print(f"완료: 결과 파일 → {self.output_file}")

    def _process_primer_pair(self, f_name: str, f_seq: str, r_name: str, r_seq: str) -> str:
        # 1) Thermo
        f_thermo = compute_thermo(f_seq)
        r_thermo = compute_thermo(r_seq)
        het_dg, het_tm = compute_heterodimer(f_seq, r_seq)

        # 2) BLAST
        blast_error = False
        f_hits_count = -1
        r_hits_count = -1
        nearby_count = -1
        min_amp_size: Optional[int] = None
        amp_details: List[str] = []

        try:
            blast_hits = run_blast_for_primers(
                f_name=f_name,
                f_seq=f_seq,
                r_name=r_name,
                r_seq=r_seq,
                blastn_path=self.blastn_path,
                db=self.blast_db,
                identity_threshold=self.identity_threshold,
                length_threshold=self.length_threshold,
                max_alignments=self.max_alignments,
            )
            f_hits_list = blast_hits.get(f_name, [])
            r_hits_list = blast_hits.get(r_name, [])
            f_hits_count = len(f_hits_list)
            r_hits_count = len(r_hits_list)

            nearby_count, min_amp_size, amp_details = find_nearby_amplicons(
                f_hits_list,
                r_hits_list,
                min_bp=self.min_amp_bp,
                max_bp=self.max_amp_bp,
            )
        except Exception:
            blast_error = True

        # 3) QC
        qc_flags, final_result, tm_diff = compute_qc_flags(
            f_thermo=f_thermo,
            r_thermo=r_thermo,
            het_dg=het_dg,
            f_hits_count=f_hits_count,
            r_hits_count=r_hits_count,
            nearby_count=nearby_count,
            blast_error=blast_error,
            tm_min=self.tm_min,
            tm_max=self.tm_max,
            tm_diff_max=self.tm_diff_max,
            hairpin_dg_cutoff=self.hairpin_dg_cutoff,
            homodimer_dg_cutoff=self.homodimer_dg_cutoff,
            heterodimer_dg_cutoff=self.heterodimer_dg_cutoff,
            blast_hit_max=self.blast_hit_max,
        )

        amplicon_info = ";".join(amp_details) if amp_details else ""

        fields = [
            f_name,
            f_seq,
            r_name,
            r_seq,
            f"{f_thermo.tm:.2f}",
            f"{r_thermo.tm:.2f}",
            f"{f_thermo.gc:.2f}",
            f"{r_thermo.gc:.2f}",
            f"{f_thermo.hairpin_dg:.2f}",
            f"{r_thermo.hairpin_dg:.2f}",
            f"{f_thermo.hairpin_tm:.2f}",
            f"{r_thermo.hairpin_tm:.2f}",
            f"{f_thermo.homodimer_dg:.2f}",
            f"{r_thermo.homodimer_dg:.2f}",
            f"{f_thermo.homodimer_tm:.2f}",
            f"{r_thermo.homodimer_tm:.2f}",
            f"{het_dg:.2f}",
            f"{het_tm:.2f}",
            str(f_hits_count),
            str(r_hits_count),
            str(nearby_count),
            "" if min_amp_size is None else str(min_amp_size),
            amplicon_info,
            qc_flags["tm_range"],
            qc_flags["tm_diff"],
            qc_flags["hairpin"],
            qc_flags["homodimer"],
            qc_flags["heterodimer"],
            qc_flags["blast_hit"],
            qc_flags["blast_amplicon"],
            final_result,
        ]
        return "\t".join(fields)

    def _write_header(self, out_handle):
        header_cols = [
            "Forward_Primer",
            "Forward_seq",
            "Reverse_Primer",
            "Reverse_seq",
            "F_Tm",
            "R_Tm",
            "F_GC%",
            "R_GC%",
            "F_Hairpin_dG_kcal",
            "R_Hairpin_dG_kcal",
            "F_Hairpin_Tm",
            "R_Hairpin_Tm",
            "F_Homodimer_dG_kcal",
            "R_Homodimer_dG_kcal",
            "F_Homodimer_Tm",
            "R_Homodimer_Tm",
            "Heterodimer_dG_kcal",
            "Heterodimer_Tm",
            "F_BLAST_hits",
            "R_BLAST_hits",
            "Nearby_amplicon_count",
            "Min_amplicon_size",
            "amplicon_info",
            "QC_Tm_range",
            "QC_Tm_diff",
            "QC_Hairpin",
            "QC_Homodimer",
            "QC_Heterodimer",
            "QC_BLAST_hit",
            "QC_BLAST_amplicon",
            "Final_Result",
        ]
        out_handle.write("\t".join(header_cols) + "\n")

def compute_thermo(seq: str) -> PrimerThermoResult:
    """단일 primer에 대해 Tm / GC / hairpin / homodimer 계산."""
    tm = primer3.calc_tm(seq)
    gc = gc_fraction(seq) * 100.0

    hairpin = primer3.calc_hairpin(seq)
    if hairpin.structure_found:
        hp_dg = hairpin.dg / 1000.0
        hp_tm = hairpin.tm
    else:
        hp_dg = 0.0
        hp_tm = 0.0

    homodimer = primer3.calc_homodimer(seq)
    if homodimer.structure_found:
        hd_dg = homodimer.dg / 1000.0
        hd_tm = homodimer.tm
    else:
        hd_dg = 0.0
        hd_tm = 0.0

    return PrimerThermoResult(
        tm=tm,
        gc=gc,
        hairpin_dg=hp_dg,
        hairpin_tm=hp_tm,
        homodimer_dg=hd_dg,
        homodimer_tm=hd_tm,
    )


def compute_heterodimer(f_seq: str, r_seq: str):
    """F/R heterodimer ΔG / Tm 계산."""
    hetero = primer3.calc_heterodimer(f_seq, r_seq)
    if hetero.structure_found:
        het_dg = hetero.dg / 1000.0
        het_tm = hetero.tm
    else:
        het_dg = 0.0
        het_tm = 0.0
    return het_dg, het_tm

# primer_qc/qc_rules.py

from typing import Dict, Tuple

from .schema import PrimerThermoResult


def compute_qc_flags(
    f_thermo: PrimerThermoResult,
    r_thermo: PrimerThermoResult,
    het_dg: float,
    f_hits_count: int,
    r_hits_count: int,
    nearby_count: int,
    blast_error: bool,
    tm_min: float,
    tm_max: float,
    tm_diff_max: float,
    hairpin_dg_cutoff: float,
    homodimer_dg_cutoff: float,
    heterodimer_dg_cutoff: float,
    blast_hit_max: int,
) -> Tuple[Dict[str, str], str, float]:
    """
    각 기준에 대해 '0' (패스) / 'x' (실패) 플래그와 최종 PASS/FAIL 반환.
    """
    # 1) Tm 범위
    qc_tm_range = (
        "0"
        if (tm_min <= f_thermo.tm <= tm_max and tm_min <= r_thermo.tm <= tm_max)
        else "x"
    )

    # 2) Tm 차이
    tm_diff = abs(f_thermo.tm - r_thermo.tm)
    qc_tm_diff = "0" if tm_diff <= tm_diff_max else "x"

    # 3) Hairpin
    qc_hairpin = (
        "0"
        if (f_thermo.hairpin_dg >= hairpin_dg_cutoff and r_thermo.hairpin_dg >= hairpin_dg_cutoff)
        else "x"
    )

    # 4) Homodimer
    qc_homodimer = (
        "0"
        if (f_thermo.homodimer_dg >= homodimer_dg_cutoff and r_thermo.homodimer_dg >= homodimer_dg_cutoff)
        else "x"
    )

    # 5) Heterodimer
    qc_heterodimer = "0" if het_dg >= heterodimer_dg_cutoff else "x"

    # 6) BLAST_hit
    if blast_error:
        qc_blast_hit = "x"
    else:
        qc_blast_hit = (
            "0" if (f_hits_count <= blast_hit_max and r_hits_count <= blast_hit_max) else "x"
        )

    # 7) BLAST_amplicon
    if blast_error:
        qc_blast_amplicon = "x"
    else:
        qc_blast_amplicon = "0" if nearby_count == 0 else "x"

    qc_flags = {
        "tm_range": qc_tm_range,
        "tm_diff": qc_tm_diff,
        "hairpin": qc_hairpin,
        "homodimer": qc_homodimer,
        "heterodimer": qc_heterodimer,
        "blast_hit": qc_blast_hit,
        "blast_amplicon": qc_blast_amplicon,
    }

    final_result = "FAIL" if "x" in qc_flags.values() else "PASS"
    return qc_flags, final_result, tm_diff