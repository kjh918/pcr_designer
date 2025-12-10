# primer_qc/blast.py

import os
import subprocess
import tempfile
from typing import Dict, List, Tuple, Optional

from .schema import PrimerBlastHit


def run_blast_for_primers(
    f_name: str,
    f_seq: str,
    r_name: str,
    r_seq: str,
    blastn_path: str,
    db: str,
    identity_threshold: float,
    length_threshold: int,
    max_alignments: int,
) -> Dict[str, List[PrimerBlastHit]]:
    """
    Forward/Reverse primer 두 개를 하나의 FASTA로 blastn-short에 넣고
    qseqid별로 hit list 반환.
    """
    fasta_lines = [f">{f_name}", f_seq, f">{r_name}", r_seq]
    fasta_str = "\n".join(fasta_lines) + "\n"

    hits: Dict[str, List[PrimerBlastHit]] = {f_name: [], r_name: []}

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "primers.fa")

        with open(fasta_path, "w") as f:
            f.write(fasta_str)

        cmd = [
            blastn_path,
            "-task", "blastn-short",
            "-db", db,
            "-query", fasta_path,
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-num_alignments",
            str(max_alignments),
        ]

        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in proc.stdout.splitlines():
            if not line.strip():
                continue

            cols = line.strip().split("\t")
            (
                qseqid,
                sseqid,
                pident,
                length,
                mismatch,
                gapopen,
                qstart,
                qend,
                sstart,
                send,
                evalue,
                bitscore,
            ) = cols

            pident = float(pident)
            length = int(length)
            qstart = int(qstart)
            qend = int(qend)
            sstart = int(sstart)
            send = int(send)
            evalue = float(evalue)
            bitscore = float(bitscore)

            # 필터 기준
            if pident < identity_threshold or length < length_threshold:
                continue

            if qseqid not in hits:
                # 예상 밖 qseqid면 스킵
                continue

            hits[qseqid].append(
                PrimerBlastHit(
                    qseqid=qseqid,
                    sseqid=sseqid,
                    pident=pident,
                    length=length,
                    qstart=qstart,
                    qend=qend,
                    sstart=sstart,
                    send=send,
                    evalue=evalue,
                    bitscore=bitscore,
                )
            )

    return hits


def hit_strand_and_3end(hit: PrimerBlastHit) -> Tuple[str, int]:
    """
    BLAST hit에서 strand(+/-)와 3' end 좌표 계산.
    sstart <= send 이면 + strand, 3' end = send
    sstart > send 이면 - strand, 3' end = sstart
    """
    if hit.sstart <= hit.send:
        strand = "+"
        three_prime = hit.send
    else:
        strand = "-"
        three_prime = hit.sstart
    return strand, three_prime


def find_nearby_amplicons(
    f_hits: List[PrimerBlastHit],
    r_hits: List[PrimerBlastHit],
    min_bp: int,
    max_bp: int,
) -> Tuple[int, Optional[int], List[str]]:
    """
    F/R hit 목록에서:
    - 같은 chr
    - strand 반대
    - 3' end 거리 min_bp~max_bp
    조합을 찾아서 (count, 최소 크기, 상세 문자열 리스트) 반환.
    """
    count = 0
    min_size: Optional[int] = None
    details: List[str] = []

    for fh in f_hits:
        f_chr = fh.sseqid
        f_strand, f_3p = hit_strand_and_3end(fh)

        for rh in r_hits:
            if rh.sseqid != f_chr:
                continue

            r_strand, r_3p = hit_strand_and_3end(rh)
            if f_strand == r_strand:
                continue

            amp_size = abs(r_3p - f_3p) + 1

            if min_bp <= amp_size <= max_bp:
                count += 1
                if min_size is None or amp_size < min_size:
                    min_size = amp_size
                details.append(f"{f_chr}:{f_3p}-{r_3p}({amp_size}bp)")

    return count, min_size, details