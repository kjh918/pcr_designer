# pcr/qc/blast.py
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


def _parse_hits(stdout: str, *, qc_params: QCParams) -> List[BlastHit]:
    """
    BLAST outfmt 6 결과를 파싱하고,
    qc_params.BLAST_IDENTITY_THRESHOLD / BLAST_LENGTH_THRESHOLD 기준으로 필터링.
    """
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
        qstart = int(cols[6])
        qend = int(cols[7])
        sstart = int(cols[8])
        send = int(cols[9])
        evalue = float(cols[10])
        bitscore = float(cols[11])
        qseq = cols[12]
        sseq = cols[13]

        # ✅ QCParams 기반 필터
        if pident < qc_params.BLAST_IDENTITY_THRESHOLD or length < qc_params.BLAST_LENGTH_THRESHOLD:
            continue

        hits.append(
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
    
def run_blast_for_single(
    name: str,
    seq: str,
    db: str,
    *,
    qc_params: QCParams,
) -> List[BlastHit]:
    """
    파일 생성 없이 표준 입력(stdin)을 이용한 단일 서열 BLAST.
    """
    name = name.strip()
    seq = seq.strip().upper()

    # 1. FASTA 형식의 문자열 생성
    fasta_str = f">{name}\n{seq}\n"

    # 2. BLAST 명령어 구성 (파일 경로 대신 '-'를 입력하여 stdin 사용 명시)
    cmd = [
        str(qc_params.BLASTN),
        "-task", "blastn-short",  # 프라이머와 같은 짧은 서열에 최적화된 옵션
        "-db", db,
        "-query", "-",            # 핵심: 표준 입력을 쿼리로 사용하도록 설정
        "-outfmt", _OUTFMT,
        "-num_alignments", str(qc_params.BLAST_MAX_ALIGNMENTS),
    ]

    # 3. subprocess 실행 (input 인자에 fasta_str을 직접 전달)
    try:
        result = subprocess.run(
            cmd, 
            input=fasta_str,      # 서열 데이터를 직접 찔러 넣음
            capture_output=True, 
            text=True, 
            check=True            # 에러 발생 시 예외 발생
        )
        return _parse_hits(result.stdout, qc_params=qc_params)
    except subprocess.CalledProcessError as e:
        # 로그 기록 등의 처리를 추가할 수 있습니다.
        print(f"BLAST 실행 에러: {e.stderr}")
        return []
    except FileNotFoundError:
        print(f"BLAST 실행 파일을 찾을 수 없습니다: {qc_params.BLASTN}")
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
    """
    primer pair를 한 번에 BLAST.
    """
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


def hit_strand_and_3end(hit: Dict[str, Any]) -> Tuple[str, int]:
    sstart = hit["sstart"]
    send = hit["send"]
    if sstart <= send:
        return "+", send
    return "-", sstart


def hit_strand_and_5end(hit: Dict[str, Any]) -> Tuple[str, int]:
    sstart = hit["sstart"]
    send = hit["send"]
    if sstart <= send:
        return "+", sstart
    return "-", send


def find_nearby_amplicons(
    f_hits: List[BlastHit],
    r_hits: List[BlastHit],
    *,
    min_bp: int,
    max_bp: int,
    f_len: Optional[int] = None,
    r_len: Optional[int] = None,
) -> Tuple[int, Optional[int], List[str]]:
    count = 0
    min_size: Optional[int] = None
    details: List[str] = []

    if f_len is None or r_len is None:
        return 0, None, []

    for fh in f_hits:
        f_chr = fh["sseqid"]
        f_strand, f_3p = hit_strand_and_3end(fh)
        _, f_5p = hit_strand_and_5end(fh)

        for rh in r_hits:
            if rh["sseqid"] != f_chr:
                continue

            r_strand, r_3p = hit_strand_and_3end(rh)

            if f_strand == r_strand:
                continue

            valid = False
            if f_strand == "+" and r_strand == "-" and f_3p < r_3p:
                valid = True
            elif f_strand == "-" and r_strand == "+" and r_3p < f_3p:
                valid = True
            if not valid:
                continue

            core_amp = abs(r_3p - f_3p) + 1
            amp_size = core_amp + f_len + r_len

            if min_bp <= amp_size <= max_bp:
                count += 1
                if min_size is None or amp_size < min_size:
                    min_size = amp_size
                details.append(
                    f"FR:{f_chr}:{f_5p}({f_strand})-{r_3p}({r_strand})({amp_size}bp)"
                )

    return count, min_size, details


def find_self_amplicons(
    hits: List[BlastHit],
    *,
    min_bp: int,
    max_bp: int,
    primer_len: Optional[int] = None,
    label: str = "F",
    primer_name: str = "PRIMER",
    primer_seq: Optional[str] = None,
) -> Tuple[int, Optional[int], List[str]]:
    count = 0
    min_size: Optional[int] = None
    details: List[str] = []

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
            if strand_i == strand_j:
                continue

            valid = False
            if strand_i == "+" and strand_j == "-" and p3_i < p3_j:
                valid = True
            elif strand_i == "-" and strand_j == "+" and p3_j < p3_i:
                valid = True
            if not valid:
                continue

            core_amp = abs(p3_j - p3_i) + 1
            amp_size = core_amp + 2 * primer_len

            if min_bp <= amp_size <= max_bp:
                count += 1
                if min_size is None or amp_size < min_size:
                    min_size = amp_size

                details.append(
                    f"{label}-SELF:{chr_i}:{p3_i}({strand_i})-{p3_j}({strand_j})({amp_size}bp)"
                    f"|pident={hi['pident']:.1f}/{hj['pident']:.1f}"
                    f"|qseq={hi['qseq']}"
                    f"|sseq_i={hi['sseq']}"
                    f"|sseq_j={hj['sseq']}"
                )

    return count, min_size, details


def probe_in_any_amplicon(
    f_hits: List[BlastHit],
    r_hits: List[BlastHit],
    p_hits: List[BlastHit],
    *,
    min_bp: int,
    max_bp: int,
    f_len: Optional[int] = None,
    r_len: Optional[int] = None,
    p_len: Optional[int] = None,
) -> Tuple[bool, List[str]]:
    if not p_hits or f_len is None or r_len is None or p_len is None:
        return False, []

    details: List[str] = []
    found = False

    for fh in f_hits:
        f_chr = fh["sseqid"]
        f_strand, f_5p = hit_strand_and_5end(fh)
        f_start = min(fh["sstart"], fh["send"])
        f_end = max(fh["sstart"], fh["send"])

        for rh in r_hits:
            if rh["sseqid"] != f_chr:
                continue

            r_strand, r_5p = hit_strand_and_5end(rh)
            r_start = min(rh["sstart"], rh["send"])
            r_end = max(rh["sstart"], rh["send"])

            if f_strand == r_strand:
                continue

            valid = False
            if f_strand == "+" and r_strand == "-" and f_5p < r_5p:
                valid = True
            elif f_strand == "-" and r_strand == "+" and r_5p < f_5p:
                valid = True
            if not valid:
                continue

            core_amp = abs(r_5p - f_5p) + 1
            amp_size = core_amp + f_len + r_len
            if not (min_bp <= amp_size <= max_bp):
                continue

            amp_start = min(f_start, r_start)
            amp_end = max(f_end, r_end)

            for ph in p_hits:
                if ph["sseqid"] != f_chr:
                    continue
                p_start = min(ph["sstart"], ph["send"])
                p_end = max(ph["sstart"], ph["send"])

                if amp_start <= p_start and p_end <= amp_end:
                    found = True
                    details.append(
                        f"PROBE_IN_FR:{f_chr}:{amp_start}-{amp_end}|probe:{p_start}-{p_end}"
                    )

    return found, details


def blast_qc_for_primer_pair(
    f_name: str,
    f_seq: str,
    r_name: str,
    r_seq: str,
    db: str,
    *,
    qc_params: QCParams,
    probe_name: Optional[str] = None,
    probe_seq: Optional[str] = None,
) -> Dict[str, Any]:
    """
    return dict 키는 기존 웹 템플릿/라우터 호환 유지:
      - f_hits, r_hits, nearby_count, min_amplicon_size, amplicon_details,
        qc_blast_hit, qc_blast_amplicon, qc_probe_in_amplicon, probe_in_amplicon, blast_error
    """
    f_seq = f_seq.strip().upper()
    r_seq = r_seq.strip().upper()
    probe_seq = probe_seq.strip().upper() if probe_seq else None

    # + 포함 시(수정염기 등) BLAST 스킵 처리(원 로직 유지)
    if "+" in f_seq or "+" in r_seq:
        return {
            "f_hits": 0,
            "r_hits": 0,
            "nearby_count": 0,
            "min_amplicon_size": None,
            "amplicon_details": [],
            "qc_blast_hit": "X",
            "qc_blast_amplicon": "X",
            "qc_probe_in_amplicon": "-" if probe_seq is None else "X",
            "probe_in_amplicon": False,
            "blast_error": False,
        }

    blast_error = False
    f_hits = r_hits = -1
    nearby_count = -1
    min_amp_size: Optional[int] = None
    amp_details: List[str] = []
    probe_in_amp = False
    probe_amp_details: List[str] = []

    min_bp = qc_params.MIN_AMP_BP
    max_bp = qc_params.MAX_AMP_BP

    try:
        blast_hits = run_blast_for_primers(f_name, f_seq, r_name, r_seq, db, qc_params=qc_params)
        f_hits_list = blast_hits.get(f_name, [])
        r_hits_list = blast_hits.get(r_name, [])

        f_hits = len(f_hits_list)
        r_hits = len(r_hits_list)

        c_FR, min_FR, det_FR = find_nearby_amplicons(
            f_hits_list,
            r_hits_list,
            min_bp=min_bp,
            max_bp=max_bp,
            f_len=len(f_seq),
            r_len=len(r_seq),
        )

        c_FF, min_FF, det_FF = find_self_amplicons(
            f_hits_list,
            min_bp=min_bp,
            max_bp=max_bp,
            primer_len=len(f_seq),
            label="F",
            primer_name=f_name,
            primer_seq=f_seq,
        )

        c_RR, min_RR, det_RR = find_self_amplicons(
            r_hits_list,
            min_bp=min_bp,
            max_bp=max_bp,
            primer_len=len(r_seq),
            label="R",
            primer_name=r_name,
            primer_seq=r_seq,
        )

        nearby_count = c_FR + c_FF + c_RR
        mins = [x for x in [min_FR, min_FF, min_RR] if x is not None]
        min_amp_size = min(mins) if mins else None
        amp_details = det_FR + det_FF + det_RR

        if probe_seq:
            pname = probe_name or "PROBE"
            p_hits_list = run_blast_for_single(pname, probe_seq, db, qc_params=qc_params)
            probe_in_amp, probe_amp_details = probe_in_any_amplicon(
                f_hits_list,
                r_hits_list,
                p_hits_list,
                min_bp=min_bp,
                max_bp=max_bp,
                f_len=len(f_seq),
                r_len=len(r_seq),
                p_len=len(probe_seq),
            )

    except Exception:
        blast_error = True

    if blast_error:
        qc_blast_hit = "X"
        qc_blast_amplicon = "X"
        qc_probe_in_amplicon = "X" if probe_seq else "-"
    else:
        # ✅ 기존 BlastQCConfig.max_hits 대체: qc_params.BLAST_MAX_ALIGNMENTS를 기준으로 동일하게 사용
        # (원래 의미가 "리포트 align 수"이긴 한데, 기존 코드의 max_hits 용도로 쓰고 있었다면 일단 동일값 사용)
        max_hits = qc_params.BLAST_MAX_ALIGNMENTS

        qc_blast_hit = "O" if (f_hits <= max_hits and r_hits <= max_hits) else "X"
        qc_blast_amplicon = "O" if nearby_count <= 1 else "X"
        qc_probe_in_amplicon = ("O" if probe_in_amp else "X") if probe_seq else "-"

    all_amp_details = amp_details + probe_amp_details

    return {
        "f_hits": f_hits,
        "r_hits": r_hits,
        "nearby_count": nearby_count,
        "min_amplicon_size": min_amp_size,
        "amplicon_details": all_amp_details,
        "qc_blast_hit": qc_blast_hit,
        "qc_blast_amplicon": qc_blast_amplicon,
        "qc_probe_in_amplicon": qc_probe_in_amplicon,
        "probe_in_amplicon": probe_in_amp,
        "blast_error": blast_error,
    }

def apply_blast_qc_to_rows(
    genomic_id: str,
    rows: List[Dict[str, Any]],
    *,
    db: str,
    qc_params: QCParams,
    include_probe: bool = True,
    fail_closed: bool = True,   # BLAST 에러 시 fail 처리 (안전)
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    rows(list[dict])에 BLAST QC 결과를 붙여서:
      - total_enriched: 전체 rows + BLAST 결과 컬럼
      - filtered: BLAST QC까지 통과한 rows

    BLAST QC PASS 조건(네 blast_qc_for_primer_pair 반환 기반):
      - qc_blast_hit == 'O'
      - qc_blast_amplicon == 'O'
      - probe 있으면 qc_probe_in_amplicon == 'O' (없으면 '-' 허용)
      - blast_error == False
    """
    total_enriched: List[Dict[str, Any]] = []
    filtered: List[Dict[str, Any]] = []

    for i, row in enumerate(rows, start=1):
        # 원본 보호 원하면 copy(); 여기선 inplace 피하려고 copy()
        a = dict(row)

        f_seq = a.get("forward_sequence") or a.get("forward_seq")
        r_seq = a.get("reverse_sequence") or a.get("reverse_seq")
        p_seq = a.get("probe_sequence") or a.get("probe_seq")

        # 이름 없으면 자동 생성 (primer pair 별로 구분되게 i 사용)
        f_name = a.get("forward_name") or f"{genomic_id}_F_{i}"
        r_name = a.get("reverse_name") or f"{genomic_id}_R_{i}"
        p_name = a.get("probe_name") or f"{genomic_id}_P_{i}"

        # 기본값(에러/결측 안전)
        a.setdefault("blast_error", False)
        a.setdefault("qc_blast_hit", "X")
        a.setdefault("qc_blast_amplicon", "X")
        a.setdefault("qc_probe_in_amplicon", "-" if not p_seq else "X")
        a.setdefault("nearby_count", -1)
        a.setdefault("min_amplicon_size", None)
        a.setdefault("f_hits", -1)
        a.setdefault("r_hits", -1)
        a.setdefault("amplicon_details", [])

        if not f_seq or not r_seq:
            # forward/reverse 없으면 BLAST 자체가 불가
            a["blast_error"] = True
            a["qc_blast_hit"] = "X"
            a["qc_blast_amplicon"] = "X"
            a["qc_probe_in_amplicon"] = "X" if (include_probe and p_seq) else "-"
        else:
            try:
                res = blast_qc_for_primer_pair(
                    f_name=f_name,
                    f_seq=str(f_seq),
                    r_name=r_name,
                    r_seq=str(r_seq),
                    db=db,
                    qc_params=qc_params,
                    probe_name=p_name if (include_probe and p_seq) else None,
                    probe_seq=str(p_seq) if (include_probe and p_seq) else None,
                )
                print(res)
                # 반환 키 그대로 merge
                a.update(res)
            except Exception:
                a["blast_error"] = True
                if fail_closed:
                    a["qc_blast_hit"] = "X"
                    a["qc_blast_amplicon"] = "X"
                    a["qc_probe_in_amplicon"] = "X" if (include_probe and p_seq) else "-"

        # ---- BLAST PASS 판정 ----
        probe_ok = (a.get("qc_probe_in_amplicon") in ("O", "-"))
        blast_pass = (
            (a.get("blast_error") is False)
            and (a.get("qc_blast_hit") == "O")
            and (a.get("qc_blast_amplicon") == "O")
            and probe_ok
        )
        a["BLAST_PASS"] = "O" if blast_pass else "X"

        total_enriched.append(a)
        if blast_pass:
            filtered.append(a)

    return total_enriched, filtered