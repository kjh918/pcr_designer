# primer/qc.py

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, List, Iterable, Any, Optional
import os
import subprocess

import primer3
from Bio.Seq import Seq

from primer.pcr_components import Primer, Amplicon
from config.settings import settings  # ⚠️ 실제 경로에 맞게 수정


# -----------------------------
# 설정에서 QC / BLAST 파라미터 로딩
# -----------------------------
_qc = settings.qc_params

# Thermo / dimer threshold
HAIRPIN_MIN_DG = _qc.HAIRPIN_MIN_DG
HOMODIMER_MIN_DG = _qc.HOMODIMER_MIN_DG
HETERODIMER_MIN_DG = _qc.HETERODIMER_MIN_DG
MAX_TM_DIFF = _qc.MAX_TM_DIFF

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

# 필요하면 qc_params에 추가해두고 쓰면 됨 (예: 허용 on-target hit 개수)
BLAST_HIT_MAX = getattr(_qc, "BLAST_HIT_MAX", 1)


# -----------------------------
# Thermo util
# -----------------------------
def compute_heterodimer(f_seq: str, r_seq: str) -> Tuple[float, float]:
    """
    두 올리고 간 heterodimer ΔG / Tm 계산.
    """
    hetero = primer3.calc_heterodimer(f_seq, r_seq)
    if hetero.structure_found:
        het_dg = hetero.dg / 1000.0  # primer3는 1000배 단위
        het_tm = hetero.tm
    else:
        het_dg = 0.0
        het_tm = 0.0
    return het_dg, het_tm


def compute_hairpin(seq: str) -> Tuple[float, float]:
    """
    단일 올리고 hairpin ΔG / Tm 계산.
    """
    hp = primer3.calc_hairpin(seq)
    if hp.structure_found:
        hp_dg = hp.dg / 1000.0
        hp_tm = hp.tm
    else:
        hp_dg = 0.0
        hp_tm = 0.0
    return hp_tm, hp_dg


def compute_homodimer(seq: str) -> float:
    """
    단일 올리고 homodimer ΔG 계산.
    """
    hd = primer3.calc_homodimer(seq)
    if hd.structure_found:
        hd_dg = hd.dg / 1000.0
    else:
        hd_dg = 0.0
    return hd_dg


# -----------------------------
# SimpleAmplicon (QC-only용)
# -----------------------------
@dataclass
class SimpleAmplicon:
    forward_sequence: str
    reverse_sequence: str
    probe_sequence: str | None = None

    def to_dict(self) -> Dict[str, Any]:
        """
        evaluate_amplicons에서 필요로 하는 최소 필드만 채운 dict 반환.
        hairpin/homodimer는 여기서 계산해서 넣어준다.
        """
        f_seq = self.forward_sequence
        r_seq = self.reverse_sequence
        p_seq = self.probe_sequence

        # Hairpin
        f_hp_tm, f_hp_dg = compute_hairpin(f_seq) if f_seq else (0.0, 0.0)
        r_hp_tm, r_hp_dg = compute_hairpin(r_seq) if r_seq else (0.0, 0.0)

        # Homodimer
        f_hd_dg = compute_homodimer(f_seq) if f_seq else 0.0
        r_hd_dg = compute_homodimer(r_seq) if r_seq else 0.0

        return {
            "forward_sequence": f_seq,
            "reverse_sequence": r_seq,
            "probe_sequence": p_seq,
            # hairpin
            "forward_hairpin_tm": f_hp_tm,
            "forward_hairpin_dg": f_hp_dg,
            "reverse_hairpin_tm": r_hp_tm,
            "reverse_hairpin_dg": r_hp_dg,
            # homodimer
            "forward_homodimer_dg": f_hd_dg,
            "reverse_homodimer_dg": r_hd_dg,
        }


# -----------------------------
# QC threshold dataclass
# -----------------------------
@dataclass
class QCThresholds:
    max_tm_diff: float = MAX_TM_DIFF
    hairpin_min_dg: float = HAIRPIN_MIN_DG
    homodimer_min_dg: float = HOMODIMER_MIN_DG
    heterodimer_min_dg: float = HETERODIMER_MIN_DG


def _qc_bool_flags(amp: Dict[str, Any], th: QCThresholds) -> Tuple[bool, bool, bool]:
    """
    내부용: hairpin / homodimer / (FR) heterodimer 각각의 True/False 리턴.
    heterodimer는 기본적으로 F–R 쌍(heterodimer_dg / heterodimer_tm)을 기준으로 판단.
    """
    # Hairpin (Tm, dG 기준)
    hairpin_ok = (
        amp.get("forward_hairpin_dg", 0.0) >= th.hairpin_min_dg
        and amp.get("reverse_hairpin_dg", 0.0) >= th.hairpin_min_dg
    )

    # Homodimer (각각 dG 기준)
    homodimer_ok = (
        amp.get("forward_homodimer_dg", 0.0) >= th.homodimer_min_dg
        and amp.get("reverse_homodimer_dg", 0.0) >= th.homodimer_min_dg
    )

    # F–R heterodimer (기본 heterodimer_dg / heterodimer_tm 사용)
    hetero_fr_ok = (
        amp.get("heterodimer_dg", 0.0) >= th.heterodimer_min_dg
        and amp.get("heterodimer_tm", 0.0) <= th.heterodimer_max_tm
    )

    return hairpin_ok, homodimer_ok, hetero_fr_ok


def amplicon_passes_qc(amp: Dict[str, Any], th: QCThresholds) -> bool:
    """
    최종 QC 통과 여부 (기본적으로 hairpin/homodimer/F–R heterodimer 기준).
    FP/RP까지 포함한 최종 판정은 evaluate_amplicons 안에서 처리.
    """
    hairpin_ok, homodimer_ok, hetero_fr_ok = _qc_bool_flags(amp, th)
    return hairpin_ok and homodimer_ok and hetero_fr_ok


def _hetero_ok(dg: float, tm: float, th: QCThresholds) -> bool:
    """
    단일 heterodimer 쌍에 대한 QC 여부.
    """
    return (dg >= th.heterodimer_min_dg) and (tm <= th.heterodimer_max_tm)


# -----------------------------
# Amplicon QC (primer3 Thermo 기반)
# -----------------------------
def evaluate_amplicons(
    genomic_id,
    amplicons: Iterable[Any],
    qc_thresholds: QCThresholds,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Amplicon 객체 리스트에 대해:
      - Amplicon.to_dict() 호출
      - F–R / F–P / R–P heterodimer_dg / heterodimer_tm 계산 추가
      - 각 QC 항목별 결과(O/X) 컬럼 및 최종 QC_PASS 컬럼 추가
      - QC 통과 여부 필터링
    """
    total_rows: List[Dict[str, Any]] = []
    filtered_rows: List[Dict[str, Any]] = []

    for amplicon in amplicons:
        a_dict = amplicon.to_dict()

        a_dict["ID"] = genomic_id
        f_seq = a_dict.get("forward_sequence")
        r_seq = a_dict.get("reverse_sequence")
        rc_r_seq = str(Seq(a_dict.get("reverse_sequence")).reverse_complement())
        a_dict["rc_reverse_sequence"] = rc_r_seq
        p_seq = a_dict.get("probe_sequence")

        # ---------- 1) heterodimer 계산 ----------
        # F–R
        if f_seq and r_seq:
            het_fr_dg, het_fr_tm = compute_heterodimer(f_seq, r_seq)
        else:
            het_fr_dg, het_fr_tm = 0.0, 0.0

        # F–P
        if f_seq and p_seq:
            het_fp_dg, het_fp_tm = compute_heterodimer(f_seq, p_seq)
        else:
            het_fp_dg, het_fp_tm = 0.0, 0.0

        # R–P
        if r_seq and p_seq:
            het_rp_dg, het_rp_tm = compute_heterodimer(r_seq, p_seq)
        else:
            het_rp_dg, het_rp_tm = 0.0, 0.0

        # 기존 heterodimer_dg / tm 은 F–R 기준으로 유지
        a_dict["heterodimer_dg"] = het_fr_dg
        a_dict["heterodimer_tm"] = het_fr_tm

        # 추가: 명시적으로 세 쌍 모두 저장
        a_dict["heterodimer_fr_dg"] = het_fr_dg
        a_dict["heterodimer_fr_tm"] = het_fr_tm
        a_dict["heterodimer_fp_dg"] = het_fp_dg
        a_dict["heterodimer_fp_tm"] = het_fp_tm
        a_dict["heterodimer_rp_dg"] = het_rp_dg
        a_dict["heterodimer_rp_tm"] = het_rp_tm

        # ---------- 2) hairpin / homodimer / (F–R) heterodimer QC ----------
        hairpin_ok, homodimer_ok, hetero_fr_ok = _qc_bool_flags(a_dict, qc_thresholds)

        # ---------- 3) F–P / R–P heterodimer QC ----------
        if f_seq and p_seq:
            hetero_fp_ok = _hetero_ok(het_fp_dg, het_fp_tm, qc_thresholds)
        else:
            hetero_fp_ok = True  # probe 없으면 이 조건은 pass 취급

        if r_seq and p_seq:
            hetero_rp_ok = _hetero_ok(het_rp_dg, het_rp_tm, qc_thresholds)
        else:
            hetero_rp_ok = True

        # ---------- 4) QC 플래그(O/X) 및 최종 QC_PASS ----------
        a_dict["QC_HAIRPIN"] = "O" if hairpin_ok else "X"
        a_dict["QC_HOMODIMER"] = "O" if homodimer_ok else "X"
        a_dict["QC_HETERODIMER_FR"] = "O" if hetero_fr_ok else "X"

        # FP/RP는 probe 유무에 따라 "O"/"X"/"-" 로 표현
        if f_seq and p_seq:
            a_dict["QC_HETERODIMER_FP"] = "O" if hetero_fp_ok else "X"
        else:
            a_dict["QC_HETERODIMER_FP"] = "-"  # probe 또는 forward 미존재

        if r_seq and p_seq:
            a_dict["QC_HETERODIMER_RP"] = "O" if hetero_rp_ok else "X"
        else:
            a_dict["QC_HETERODIMER_RP"] = "-"  # probe 또는 reverse 미존재

        qc_pass = (
            hairpin_ok
            and homodimer_ok
            and hetero_fr_ok
            and hetero_fp_ok
            and hetero_rp_ok
        )
        a_dict["QC_PASS"] = "O" if qc_pass else "X"

        total_rows.append(a_dict)
        if qc_pass:
            filtered_rows.append(a_dict)

    return total_rows, filtered_rows


# -----------------------------
# BLAST 기반 QC
# -----------------------------
@dataclass
class BlastQCConfig:
    identity_threshold: float = BLAST_IDENTITY_THRESHOLD
    length_threshold: int = BLAST_LENGTH_THRESHOLD
    min_amp_bp: int = MIN_AMP_BP
    max_amp_bp: int = MAX_AMP_BP
    max_hits: int = BLAST_HIT_MAX
    max_alignments: int = BLAST_MAX_ALIGNMENTS


def hit_strand_and_3end(hit: Dict[str, Any]) -> Tuple[str, int]:
    """
    BLAST hit에서 strand(+/-)와 3' end 좌표를 계산.
    """
    sstart = hit["sstart"]
    send = hit["send"]
    if sstart <= send:
        strand = "+"
        three_prime = send
    else:
        strand = "-"
        three_prime = sstart
    return strand, three_prime


def find_nearby_amplicons(
    f_hits,
    r_hits,
    min_bp=MIN_AMP_BP,
    max_bp=MAX_AMP_BP,
    f_len=None,
    r_len=None,
):
    """
    서로 다른 primer(F/R) 사이의 잠재 amplicon.
    PCR product 가능 조합만 카운트:
    - 같은 chr
    - 서로 반대 strand
    - 서로를 향하는 방향
    - 길이: core distance + F/R primer 길이
    """
    count = 0
    min_size = None
    details = []

    if f_len is None or r_len is None:
        return 0, None, []

    for fh in f_hits:
        f_chr = fh["sseqid"]
        f_strand, f_3p = hit_strand_and_3end(fh)

        for rh in r_hits:
            if rh["sseqid"] != f_chr:
                continue

            r_strand, r_3p = hit_strand_and_3end(rh)

            # 서로 다른 strand
            if f_strand == r_strand:
                continue

            # 서로를 향하는 방향인지 체크
            valid_orientation = False
            if f_strand == "+" and r_strand == "-" and f_3p < r_3p:
                valid_orientation = True
            elif f_strand == "-" and r_strand == "+" and r_3p < f_3p:
                valid_orientation = True

            if not valid_orientation:
                continue

            core_amp = abs(r_3p - f_3p) + 1
            amp_size = core_amp + f_len + r_len

            if min_bp <= amp_size <= max_bp:
                count += 1
                if min_size is None or amp_size < min_size:
                    min_size = amp_size
                details.append(
                    f"FR:{f_chr}:{f_3p}({f_strand})-{r_3p}({r_strand})({amp_size}bp)"
                )

    return count, min_size, details


def find_self_amplicons(
    hits,
    min_bp=MIN_AMP_BP,
    max_bp=MAX_AMP_BP,
    primer_len=None,
    label="F",
    primer_name="PRIMER",
    primer_seq=None,
    db=None,
    plot_dir=None,
):
    """
    같은 primer(hits) 안에서 self-amplicon 찾기
    - 같은 chr
    - 서로 반대 strand
    - 서로를 향하는 방향
    - 길이: core distance + 2 * primer_len
    """
    count = 0
    min_size = None
    details = []

    if primer_len is None or primer_seq is None:
        return 0, None, []

    n = len(hits)
    for i in range(n):
        hi = hits[i]
        chr_i = hi["sseqid"]
        strand_i, p3_i = hit_strand_and_3end(hi)

        for j in range(i + 1, n):
            hj = hits[j]
            if hj["sseqid"] != chr_i:
                continue

            strand_j, p3_j = hit_strand_and_3end(hj)

            # 서로 다른 strand
            if strand_i == strand_j:
                continue

            # 서로를 향하는 방향인지 체크
            valid_orientation = False
            if strand_i == "+" and strand_j == "-" and p3_i < p3_j:
                valid_orientation = True
            elif strand_i == "-" and strand_j == "+" and p3_j < p3_i:
                valid_orientation = True

            if not valid_orientation:
                continue

            core_amp = abs(p3_j - p3_i) + 1
            amp_size = core_amp + 2 * primer_len

            if min_bp <= amp_size <= max_bp:
                count += 1
                if min_size is None or amp_size < min_size:
                    min_size = amp_size

                pident_i = hi["pident"]
                pident_j = hj["pident"]
                qseq = hi["qseq"]
                sseq_i = hi["sseq"]
                sseq_j = hj["sseq"]

                info = (
                    f"{label}-SELF:{chr_i}:{p3_i}({strand_i})-{p3_j}({strand_j})({amp_size}bp)"
                    f"|pident={pident_i:.1f}/{pident_j:.1f}"
                    f"|qseq={qseq}"
                    f"|sseq_i={sseq_i}"
                    f"|sseq_j={sseq_j}"
                )
                details.append(info)

                # QC-only에서는 plot_dir=None으로 넘기면 이 부분은 실행 안 됨
                # if db is not None and plot_dir is not None:
                #     ...

    return count, min_size, details


def run_blast_for_primers(f_name, f_seq, r_name, r_seq, db):
    """
    Forward, Reverse 프라이머를 하나의 FASTA로 만들어 blastn 수행 후
    각 primer별 hit 목록을 반환.
    hit dict에 qseq/sseq까지 포함.
    """
    orig_f_name = f_name  # BLAST qseqid용
    orig_r_name = r_name

    safe_f_name = (
        orig_f_name.replace("(", "_")
        .replace(")", "_")
        .replace(" ", "_")
        .replace("/", "_")
    )
    fasta_str = f">{orig_f_name}\n{f_seq}\n>{orig_r_name}\n{r_seq}\n"

    fasta_path = f"{safe_f_name}.fasta"
    with open(fasta_path, "w") as tmp_fasta:
        tmp_fasta.write(fasta_str)

    outfmt = (
        "6 qseqid sseqid pident length mismatch gapopen qstart qend "
        "sstart send evalue bitscore qseq sseq"
    )

    cmd = [
        BLASTN,
        "-task",
        "blastn-short",
        "-db",
        db,
        "-query",
        fasta_path,
        "-outfmt",
        outfmt,
        "-num_alignments",
        str(BLAST_MAX_ALIGNMENTS),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return {orig_f_name: [], orig_r_name: []}

    hits = {
        orig_f_name: [],
        orig_r_name: [],
    }

    for line in result.stdout.strip().splitlines():
        if not line.strip():
            continue
        cols = line.strip("\n").split("\t")
        if len(cols) < 14:
            continue

        qseqid = cols[0]
        sseqid = cols[1]
        pident = float(cols[2])
        length = int(cols[3])
        qstart = int(cols[6])
        qend = int(cols[7])
        sstart = int(cols[8])
        send = int(cols[9])
        evalue = float(cols[10])
        bitscore = float(cols[11])
        qseq = cols[12]
        sseq = cols[13]

        # 필터 기준 적용 (identity/length는 config에도 있지만,
        # 여기서는 상단에서 settings로부터 불러온 상수를 사용)
        if pident < BLAST_IDENTITY_THRESHOLD or length < BLAST_LENGTH_THRESHOLD:
            continue

        if qseqid not in hits:
            continue

        hits[qseqid].append(
            {
                "qseqid": qseqid,
                "sseqid": sseqid,
                "pident": pident,
                "length": length,
                "qstart": qstart,
                "qend": qend,
                "sstart": sstart,
                "send": send,
                "evalue": evalue,
                "bitscore": bitscore,
                "qseq": qseq,
                "sseq": sseq,
            }
        )

    return hits


def blast_qc_for_primer_pair(
    f_name: str,
    f_seq: str,
    r_name: str,
    r_seq: str,
    db: str,
    *,
    config: BlastQCConfig = BlastQCConfig(),
) -> Dict[str, Any]:
    """
    Forward / Reverse primer 한 쌍에 대해 BLAST 기반 QC 계산.
    """
    f_seq = f_seq.strip().upper()
    r_seq = r_seq.strip().upper()

    # '+' 포함된 primer는 바로 FAIL 또는 스킵 처리
    if "+" in f_seq or "+" in r_seq:
        return {
            "f_hits": 0,
            "r_hits": 0,
            "nearby_count": 0,
            "min_amplicon_size": None,
            "amplicon_details": [],
            "qc_blast_hit": "X",
            "qc_blast_amplicon": "X",
            "blast_error": False,
        }

    blast_error = False
    f_hits = r_hits = -1
    nearby_count = -1
    min_amp_size: Optional[int] = None
    amp_details: List[str] = []

    try:
        blast_hits = run_blast_for_primers(f_name, f_seq, r_name, r_seq, db)
        f_hits_list = blast_hits.get(f_name, [])
        r_hits_list = blast_hits.get(r_name, [])

        f_hits = len(f_hits_list)
        r_hits = len(r_hits_list)

        # F-R amplicon
        c_FR, min_FR, det_FR = find_nearby_amplicons(
            f_hits_list,
            r_hits_list,
            min_bp=config.min_amp_bp,
            max_bp=config.max_amp_bp,
            f_len=len(f_seq),
            r_len=len(r_seq),
        )

        # F-SELF
        c_FF, min_FF, det_FF = find_self_amplicons(
            f_hits_list,
            min_bp=config.min_amp_bp,
            max_bp=config.max_amp_bp,
            primer_len=len(f_seq),
            label="F",
            primer_name=f_name,
            primer_seq=f_seq,
            db=db,
            plot_dir=None,  # QC-only에서는 그림 안 그림
        )

        # R-SELF
        c_RR, min_RR, det_RR = find_self_amplicons(
            r_hits_list,
            min_bp=config.min_amp_bp,
            max_bp=config.max_amp_bp,
            primer_len=len(r_seq),
            label="R",
            primer_name=r_name,
            primer_seq=r_seq,
            db=db,
            plot_dir=None,
        )

        nearby_count = c_FR + c_FF + c_RR
        mins = [x for x in [min_FR, min_FF, min_RR] if x is not None]
        min_amp_size = min(mins) if mins else None
        amp_details = det_FR + det_FF + det_RR

    except Exception:
        blast_error = True

    if blast_error:
        qc_blast_hit = "X"
        qc_blast_amplicon = "X"
    else:
        qc_blast_hit = (
            "O" if (f_hits <= config.max_hits and r_hits <= config.max_hits) else "X"
        )
        qc_blast_amplicon = "O" if nearby_count == 0 else "X"

    return {
        "f_hits": f_hits,
        "r_hits": r_hits,
        "nearby_count": nearby_count,
        "min_amplicon_size": min_amp_size,
        "amplicon_details": amp_details,
        "qc_blast_hit": qc_blast_hit,
        "qc_blast_amplicon": qc_blast_amplicon,
        "blast_error": blast_error,
    }


# -----------------------------
# QC-only용 Amplicon 생성
# -----------------------------
def make_amplicon_for_qc(
    forward_seq: str,
    reverse_seq: str,
    probe_seq: str | None = None,
    template_seq: str | None = None,
) -> Amplicon:
    """
    QC only용 Amplicon 생성:
    - template_sequence가 주어지면 그걸 사용
    - 없으면 forward + 'N'*10 + reverse 를 이어붙인 가짜 템플릿 생성
    """
    forward_seq = forward_seq.strip().upper()
    reverse_seq = reverse_seq.strip().upper()
    probe_seq = probe_seq.strip().upper() if probe_seq else None

    if template_seq and template_seq.strip():
        template = template_seq.strip().upper()
    else:
        template = forward_seq + ("N" * 10) + reverse_seq

    target_start_index = 0
    target_end_index = len(template) - 1

    f_primer = Primer(
        template_sequence=template,
        reference_template_sequence=template,
        sequence=forward_seq,
        strand="forward",
        primer_type="forward",
        target_start_index=target_start_index,
        target_end_index=target_end_index,
    )

    r_primer = Primer(
        template_sequence=template,
        reference_template_sequence=template,
        sequence=reverse_seq,
        strand="reverse",
        primer_type="reverse",
        target_start_index=target_start_index,
        target_end_index=target_end_index,
    )

    probe_primer: Primer | None = None
    if probe_seq:
        probe_primer = Primer(
            template_sequence=template,
            reference_template_sequence=template,
            sequence=probe_seq,
            strand="forward",
            primer_type="probe",
            target_start_index=target_start_index,
            target_end_index=target_end_index,
        )

    amp = Amplicon(
        template_sequence=template,
        reference_template_sequence=template,
        target_start_index=target_start_index,
        target_end_index=target_end_index,
        forward_primer=f_primer,
        reverse_primer=r_primer,
        probe=probe_primer,
    )

    return amp
