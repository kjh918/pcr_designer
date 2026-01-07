from __future__ import annotations

from typing import Dict, Any, List, Tuple, Optional
import subprocess
import tempfile
import os

from pcr.config.schema.qc import QCParams
from pcr.qc.types import BlastHit


_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen qstart qend "
    "sstart send evalue bitscore qseq sseq"
)


# ---------------------------
# 0) 작은 유틸
# ---------------------------
def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def hit_strand(hit: BlastHit) -> str:
    return "+" if hit["sstart"] <= hit["send"] else "-"


def hit_subject_interval_1b(hit: BlastHit) -> Tuple[int, int]:
    s = min(hit["sstart"], hit["send"])
    e = max(hit["sstart"], hit["send"])
    return s, e


def hit_3p_1b(hit: BlastHit) -> int:
    s, e = hit_subject_interval_1b(hit)
    return e if hit_strand(hit) == "+" else s


def _fetch_seq_from_fasta(
    fasta_path: Optional[str],
    chrom: str,
    start0: int,
    end0_excl: int,
) -> Optional[str]:
    if not fasta_path:
        return None
    try:
        import pysam  # type: ignore
    except Exception:
        return None

    start0 = max(0, int(start0))
    end0_excl = max(start0, int(end0_excl))
    try:
        fa = pysam.FastaFile(str(fasta_path))
        seq = fa.fetch(chrom, start0, end0_excl)
        fa.close()
        seq = (seq or "").upper()
        return seq if seq else None
    except Exception:
        return None


def _parse_hits(stdout: str, *, qc_params: QCParams) -> List[BlastHit]:
    hits: List[BlastHit] = []
    for line in stdout.strip().splitlines():
        if not line.strip():
            continue
        cols = line.strip("\n").split("\t")
        if len(cols) < 14:
            continue

        qseqid = cols[0]
        sseqid = cols[1]
        pident = float(cols[2])
        length = int(cols[3])
        mismatch = int(cols[4])
        gapopen = int(cols[5])
        qstart = int(cols[6])
        qend = int(cols[7])
        sstart = int(cols[8])
        send = int(cols[9])
        evalue = float(cols[10])
        bitscore = float(cols[11])
        qseq = cols[12]
        sseq = cols[13]

        # ✅ QC 통과(필터링) 조건
        if pident < qc_params.BLAST_IDENTITY_THRESHOLD or length < qc_params.BLAST_LENGTH_THRESHOLD:
            continue

        hits.append(
            {
                "qseqid": qseqid,
                "sseqid": sseqid,
                "pident": pident,
                "length": length,
                "mismatch": mismatch,
                "gapopen": gapopen,
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


def run_blast_for_single(name: str, seq: str, db: str, *, qc_params: QCParams) -> List[BlastHit]:
    name = name.strip()
    seq = seq.strip().upper()
    fasta_str = f">{name}\n{seq}\n"

    cmd = [
        str(qc_params.BLASTN),
        "-task", "blastn-short",
        "-db", db,
        "-query", "-",
        "-outfmt", _OUTFMT,
        "-num_alignments", str(qc_params.BLAST_MAX_ALIGNMENTS),
    ]

    try:
        result = subprocess.run(
            cmd,
            input=fasta_str,
            capture_output=True,
            text=True,
            check=True,
        )
        return _parse_hits(result.stdout, qc_params=qc_params)
    except Exception:
        return []


def run_blast_for_primers(
    f_name: str,
    f_seq: str,
    r_name: str,
    r_seq: str,
    db: str,
    *,
    qc_params: QCParams,
) -> Dict[str, List[BlastHit]]:
    f_seq = f_seq.strip().upper()
    r_seq = r_seq.strip().upper()
    fasta_str = f">{f_name}\n{f_seq}\n>{r_name}\n{r_seq}\n"

    with tempfile.TemporaryDirectory() as td:
        fasta_path = os.path.join(td, "query.fasta")
        with open(fasta_path, "w") as f:
            f.write(fasta_str)

        cmd = [
            str(qc_params.BLASTN),
            "-task", "blastn-short",
            "-db", db,
            "-query", fasta_path,
            "-outfmt", _OUTFMT,
            "-num_alignments", str(qc_params.BLAST_MAX_ALIGNMENTS),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return {f_name: [], r_name: []}

        all_hits = _parse_hits(result.stdout, qc_params=qc_params)
        hits: Dict[str, List[BlastHit]] = {f_name: [], r_name: []}
        for h in all_hits:
            if h["qseqid"] in hits:
                hits[h["qseqid"]].append(h)
        return hits


# ---------------------------
# 1) QC 통과된 hit들로 FR 페어(amplicon 생성가능) 찾기
# ---------------------------
def _is_valid_fr_pair(fh: BlastHit, rh: BlastHit) -> bool:
    if fh["sseqid"] != rh["sseqid"]:
        return False
    fs = hit_strand(fh)
    rs = hit_strand(rh)
    if fs == rs:
        return False

    f3 = hit_3p_1b(fh)
    r3 = hit_3p_1b(rh)

    if fs == "+" and rs == "-" and f3 < r3:
        return True
    if fs == "-" and rs == "+" and r3 < f3:
        return True
    return False


# ✅ ADDED: 여러 FR 페어를 "ref window 길이" 기준으로 정렬해서 반환
def _build_all_fr_pairs(
    f_hits: List[BlastHit],
    r_hits: List[BlastHit],
    *,
    min_bp: int,
    max_bp: int,
    f_len: int,
    r_len: int,
    max_pairs: int = 50,
) -> List[Tuple[BlastHit, BlastHit, int]]:
    """
    return: [(fh, rh, ref_len), ...] sorted by ref_len asc
    """
    pairs: List[Tuple[BlastHit, BlastHit, int]] = []

    for fh in f_hits:
        for rh in r_hits:
            if not _is_valid_fr_pair(fh, rh):
                continue

            f3 = hit_3p_1b(fh)
            r3 = hit_3p_1b(rh)
            fs = hit_strand(fh)
            rs = hit_strand(rh)

            # forward full span (1-based inclusive)
            if fs == "+":
                f_full_s = f3 - f_len + 1
                f_full_e = f3
            else:
                f_full_s = f3
                f_full_e = f3 + f_len - 1

            # reverse full span
            if rs == "+":
                r_full_s = r3 - r_len + 1
                r_full_e = r3
            else:
                r_full_s = r3
                r_full_e = r3 + r_len - 1

            ref_s = min(f_full_s, r_full_s)
            ref_e = max(f_full_e, r_full_e)
            ref_len = ref_e - ref_s + 1

            if not (min_bp <= ref_len <= max_bp):
                continue

            pairs.append((fh, rh, ref_len))

    pairs.sort(key=lambda x: x[2])
    return pairs[:max_pairs]


# ---------------------------
# 2) reference sequence 만들기 (primer 길이 고려)
# ---------------------------
def _build_reference_window(
    *,
    chrom: str,
    fh: BlastHit,
    rh: BlastHit,
    f_seq: str,
    r_seq: str,
) -> Dict[str, Any]:
    f_len = len(f_seq)
    r_len = len(r_seq)

    f3 = hit_3p_1b(fh)
    r3 = hit_3p_1b(rh)
    fs = hit_strand(fh)
    rs = hit_strand(rh)

    if fs == "+":
        f_full_s = f3 - f_len + 1
        f_full_e = f3
    else:
        f_full_s = f3
        f_full_e = f3 + f_len - 1

    if rs == "+":
        r_full_s = r3 - r_len + 1
        r_full_e = r3
    else:
        r_full_s = r3
        r_full_e = r3 + r_len - 1

    ref_s = min(f_full_s, r_full_s)
    ref_e = max(f_full_e, r_full_e)

    return {
        "chrom": chrom,
        "reference_start_1b": ref_s,
        "reference_end_1b": ref_e,
        "forward_full_bind_1b": (f_full_s, f_full_e),
        "reverse_full_bind_1b": (r_full_s, r_full_e),
    }


# ---------------------------
# 3) reference 기준으로 amplicon(사이 영역) + mismatch index 만들기
# ---------------------------
def _compare_primer_to_reference(
    *,
    reference_seq: str,
    reference_start_1b: int,
    primer_seq: str,
    primer_strand: str,
    full_bind_1b: Tuple[int, int],
    hit: BlastHit,
) -> Dict[str, Any]:
    primer_seq = primer_seq.upper()
    L = len(primer_seq)

    b_s, b_e = full_bind_1b
    s0 = b_s - reference_start_1b
    e0_excl = (b_e - reference_start_1b) + 1
    ref_bind = reference_seq[s0:e0_excl]

    ref_for_compare = _revcomp(ref_bind) if primer_strand == "-" else ref_bind

    mm_primer_idx: List[int] = []
    mm_detail: List[Dict[str, Any]] = []

    for i in range(min(L, len(ref_for_compare))):
        if primer_seq[i] != ref_for_compare[i]:
            mm_primer_idx.append(i)
            ref_i0 = (s0 + i) if primer_strand == "+" else ((e0_excl - 1) - i)
            mm_detail.append({
                "primer_i0": i,
                "ref_i0": ref_i0,
                "ref_pos_1b": reference_start_1b + ref_i0,
                "primer_base": primer_seq[i],
                "ref_base": reference_seq[ref_i0],
            })

    qstart = int(hit.get("qstart", 1))
    qend = int(hit.get("qend", 0))
    qstart = max(1, min(L, qstart))
    qend = max(0, min(L, qend))

    unmatched_5p = list(range(0, max(0, qstart - 1)))
    unmatched_3p = list(range(min(L, qend), L))  # ✅ 말단

    mismatch_all = sorted(set(mm_primer_idx + unmatched_5p + unmatched_3p))

    return {
        "primer_seq": primer_seq,
        "strand": primer_strand,
        "full_bind_1b": {"start": b_s, "end": b_e},
        "ref_bind_seq": ref_bind,
        "qstart": qstart,
        "qend": qend,
        "unmatched_5p_indices": unmatched_5p,
        "unmatched_3p_indices": unmatched_3p,
        "mismatch_aligned_indices": sorted(set(mm_primer_idx)),
        "mismatch_all_indices": mismatch_all,
        "mismatch_detail": mm_detail,
    }


# ✅ ADDED: 단일 (fh,rh)로부터 "amplicon 후보 1개" 결과 dict 생성
def _build_one_amplicon_result(
    *,
    fh: BlastHit,
    rh: BlastHit,
    f_seq: str,
    r_seq: str,
    fasta_path: str,
    qc_params: QCParams,
    db: str,
    probe_name: Optional[str],
    probe_seq: Optional[str],
) -> Optional[Dict[str, Any]]:
    chrom = fh["sseqid"]

    ref_meta = _build_reference_window(
        chrom=chrom,
        fh=fh,
        rh=rh,
        f_seq=f_seq,
        r_seq=r_seq,
    )
    ref_s = int(ref_meta["reference_start_1b"])
    ref_e = int(ref_meta["reference_end_1b"])

    reference_seq = _fetch_seq_from_fasta(fasta_path, chrom, ref_s - 1, ref_e)
    if not reference_seq:
        return None

    # inner amplicon (3'~3')
    f3 = hit_3p_1b(fh)
    r3 = hit_3p_1b(rh)

    if hit_strand(fh) == "+" and hit_strand(rh) == "-":
        amp_s = f3 + 1
        amp_e = r3 - 1
    else:
        amp_s = r3 + 1
        amp_e = f3 - 1

    if amp_s > amp_e:
        amplicon_seq = ""
        amplicon_len = 0
    else:
        s0 = amp_s - ref_s
        e0_excl = (amp_e - ref_s) + 1
        amplicon_seq = reference_seq[s0:e0_excl]
        amplicon_len = len(amplicon_seq)

    f_cmp = _compare_primer_to_reference(
        reference_seq=reference_seq,
        reference_start_1b=ref_s,
        primer_seq=f_seq,
        primer_strand=hit_strand(fh),
        full_bind_1b=tuple(ref_meta["forward_full_bind_1b"]),
        hit=fh,
    )
    r_cmp = _compare_primer_to_reference(
        reference_seq=reference_seq,
        reference_start_1b=ref_s,
        primer_seq=r_seq,
        primer_strand=hit_strand(rh),
        full_bind_1b=tuple(ref_meta["reverse_full_bind_1b"]),
        hit=rh,
    )

    # probe in inner amplicon?
    probe_in_amp = None
    if probe_seq:
        probe_in_amp = False
        p_hits = run_blast_for_single(probe_name or "PROBE", probe_seq, db, qc_params=qc_params)
        for ph in p_hits:
            if ph["sseqid"] != chrom:
                continue
            ps, pe = hit_subject_interval_1b(ph)
            if amp_s <= ps and pe <= amp_e:
                probe_in_amp = True
                break

    # PASS: probe 조건만 기본 반영(원하면 mismatch 조건도 여기 추가)
    PASS = True
    if probe_seq and probe_in_amp is not True:
        PASS = False

    return {
        "PASS": PASS,
        "genomic": {
            "chrom": chrom,
            "forward_hit_interval_1b": {"start": hit_subject_interval_1b(fh)[0], "end": hit_subject_interval_1b(fh)[1]},
            "reverse_hit_interval_1b": {"start": hit_subject_interval_1b(rh)[0], "end": hit_subject_interval_1b(rh)[1]},
            "forward_strand": hit_strand(fh),
            "reverse_strand": hit_strand(rh),
        },
        "reference": {
            "start_1b": ref_s,
            "end_1b": ref_e,
            "seq": reference_seq,
            "forward_full_bind_1b": {"start": ref_meta["forward_full_bind_1b"][0], "end": ref_meta["forward_full_bind_1b"][1]},
            "reverse_full_bind_1b": {"start": ref_meta["reverse_full_bind_1b"][0], "end": ref_meta["reverse_full_bind_1b"][1]},
        },
        "amplicon": {
            "start_1b": amp_s,
            "end_1b": amp_e,
            "length_bp": amplicon_len,
            "seq": amplicon_seq,
        },
        "binding": {
            "forward": f_cmp,
            "reverse": r_cmp,
            "probe_in_amplicon": probe_in_amp if probe_seq else None,
        },
    }


def blast_qc_for_primer_pair(
    f_name: str,
    f_seq: str,
    r_name: str,
    r_seq: str,
    db: str,
    *,
    fasta: str,
    qc_params: QCParams,
    probe_name: Optional[str] = None,
    probe_seq: Optional[str] = None,
    max_amplicons: int = 20,  # ✅ ADDED: 생성 후보 제한
) -> Dict[str, Any]:
    """
    ✅ CHANGED: 여러 amplicon 후보를 만들어서 반환
      - result["amplicons"] = [amplicon_result1, amplicon_result2, ...]
      - PASS = 후보 중 PASS=True 가 하나라도 있으면 True
      - filtered_amplicons = PASS=True 인 것만
    """
    f_seq = f_seq.strip().upper()
    r_seq = r_seq.strip().upper()
    probe_seq = probe_seq.strip().upper() if probe_seq else None

    if "+" in f_seq or "+" in r_seq:
        return {"blast_error": False, "PASS": False, "result": {"amplicons": [], "filtered_amplicons": []}}

    fasta_path = str(fasta[0]) if isinstance(fasta, (tuple, list)) else str(fasta)

    try:
        blast_hits = run_blast_for_primers(f_name, f_seq, r_name, r_seq, db, qc_params=qc_params)
        f_hits = blast_hits.get(f_name, [])
        r_hits = blast_hits.get(r_name, [])

        # ✅ ADDED: 가능한 FR 페어들을 모두 구성
        pairs = _build_all_fr_pairs(
            f_hits,
            r_hits,
            min_bp=qc_params.MIN_AMP_BP,
            max_bp=qc_params.MAX_AMP_BP,
            f_len=len(f_seq),
            r_len=len(r_seq),
            max_pairs=max_amplicons,
        )

        amplicons: List[Dict[str, Any]] = []
        for fh, rh, _ref_len in pairs:
            one = _build_one_amplicon_result(
                fh=fh,
                rh=rh,
                f_seq=f_seq,
                r_seq=r_seq,
                fasta_path=fasta_path,
                qc_params=qc_params,
                db=db,
                probe_name=probe_name,
                probe_seq=probe_seq,
            )
            if one is not None:
                amplicons.append(one)

        filtered = [a for a in amplicons if a.get("PASS") is True]
        PASS = len(filtered) > 0

        return {
            "blast_error": False,
            "PASS": PASS,
            "result": {
                "amplicons": amplicons,
                "filtered_amplicons": filtered,  # ✅ QC 진행은 이걸로 하면 됨
            },
        }

    except Exception as e:
        print(f"[blast_qc_for_primer_pair] error: {e}")
        return {"blast_error": True, "PASS": False, "result": {"amplicons": [], "filtered_amplicons": []}}
