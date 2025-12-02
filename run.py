import os
import subprocess
import tempfile

import primer3
from Bio.SeqUtils import gc_fraction

# --------- 설정 부분 ---------
input_path = '/storage/home/jhkim/Projects/Task/GCX-HearlingLoss_PrimerValidation-2025-12-02/Resources/primer_thermo_result.remove_gc_clamp.tsv'
output_dir = '/storage/home/jhkim/Projects/Task/GCX-HearlingLoss_PrimerValidation-2025-12-02/Results'
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'primer_thermo_blast_result.tab')

# BLAST DB (makeblastdb로 만든 DB 이름/경로)
blast_db = '/storage/home/jhkim/Apps/ncbi-blast-2.17.0+/ref/hg38'  # TODO: 여기 실제 DB 경로로 수정

# BLAST 필터 기준
BLAST_IDENTITY_THRESHOLD = 85.0   # % (조금 여유 있게 80으로)
BLAST_LENGTH_THRESHOLD = 12       # 최소 매칭 길이 (nt)
BLAST_MAX_ALIGNMENTS = 200        # 최대 리포트 align 수

# Thermo 기준 (kcal/mol)
HAIRPIN_DG_CUTOFF = -5.0
HOMODIMER_DG_CUTOFF = -6.0
HETERODIMER_DG_CUTOFF = -6.0

# Tm 기준
TM_MIN = 55.0
TM_MAX = 70.0
TM_DIFF_MAX = 3.0

# BLAST hit 개수 기준
BLAST_HIT_MAX = 1  # 자기 타겟 1개만 있다고 보는 보수적 기준

# F/R 같은 chr 내 잠재 amplicon 거리 기준 (bp)
MIN_AMP_BP = 50
MAX_AMP_BP = 300


def run_blast_for_primers(f_name, f_seq, r_name, r_seq, db):
    """
    Forward, Reverse 프라이머를 하나의 FASTA로 만들어 blastn 수행 후
    각 primer별 hit 목록을 반환.
    hit는 dict: {
        'qseqid', 'sseqid', 'pident', 'length',
        'qstart', 'qend', 'sstart', 'send',
        'evalue', 'bitscore'
    }
    """
    fasta_lines = [
        f">{f_name}",
        f_seq,
        f">{r_name}",
        r_seq
    ]
    fasta_str = "\n".join(fasta_lines) + "\n"

    with open('temp.fasta','w') as tmp_fasta:
        tmp_fasta.write(fasta_str)
    cmd = [
        "/storage/home/jhkim/Apps/ncbi-blast-2.17.0+/bin/blastn",
        "-task", "blastn-short",
        "-db", db,
        "-query", 'temp.fasta',
        # qseqid sseqid까지 포함해서 받기
        "-outfmt", "6",
        "-num_alignments", str(BLAST_MAX_ALIGNMENTS),
        ">", 'temp.txt',
    ]
    os.system(' '.join(cmd))

    hits = {
        f_name: [],
        r_name: [],
    }
    with open('temp.txt', 'r') as splitlines:
        for line in splitlines:
            cols = line.strip('\n').split("\t")
            qseqid = cols[0]
            sseqid = cols[1]
            pident = float(cols[2])
            length = int(cols[3])
            # mismatch = int(cols[4])
            # gapopen = int(cols[5])
            qstart = int(cols[6])
            qend = int(cols[7])
            sstart = int(cols[8])
            send = int(cols[9])
            evalue = float(cols[10])
            bitscore = float(cols[11])

            # 필터 기준 적용
            if pident < BLAST_IDENTITY_THRESHOLD or length < BLAST_LENGTH_THRESHOLD:
                continue

            if qseqid not in hits:
                # 예상 밖 qseqid면 스킵
                continue

            hits[qseqid].append({
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
            })

    return hits


def hit_strand_and_3end(hit):
    """
    BLAST hit에서 strand(+/-)와 3' end 좌표를 계산.
    sstart <= send 이면 + strand, 3' end = send
    sstart > send 이면 - strand, 3' end = sstart
    """
    sstart = hit["sstart"]
    send = hit["send"]
    if sstart <= send:
        strand = '+'
        three_prime = send
    else:
        strand = '-'
        three_prime = sstart
    return strand, three_prime

def find_nearby_amplicons(f_hits, r_hits,
                          min_bp=MIN_AMP_BP,
                          max_bp=MAX_AMP_BP):
    """
    F/R hit 목록에서,
    - 같은 sseqid(chr)이고
    - 서로 반대 strand이고
    - 3' end 거리(amp size)가 min_bp~max_bp 사이인 조합을 찾음.

    반환:
    (count, min_size, details_list)
    details 예: "chr21:42382147(+)-42382127(-)(123bp)"
    """
    count = 0
    min_size = None
    details = []

    for fh in f_hits:
        f_chr = fh["sseqid"]
        f_strand, f_3p = hit_strand_and_3end(fh)

        for rh in r_hits:
            if rh["sseqid"] != f_chr:
                continue

            r_strand, r_3p = hit_strand_and_3end(rh)

            # 보통 PCR은 F/R 반대 strand여야 함
            if f_strand == r_strand:
                continue

            amp_size = abs(r_3p - f_3p) + 1

            if min_bp <= amp_size <= max_bp:
                count += 1
                if min_size is None or amp_size < min_size:
                    min_size = amp_size
                details.append(
                    f"{f_chr}:{f_3p}-{r_3p}({amp_size}bp)"
                )

    return count, min_size, details


with open(input_path, 'r') as handle, open(output_file, 'w') as out:

    # 결과 헤더 (QC 컬럼은 0/ x, amplicon_info 추가)
    out.write(
        '\t'.join([
            'Forward_Primer', 'Forward_seq',
            'Reverse_Primer', 'Reverse_seq',
            'F_Tm', 'R_Tm',
            'F_GC%', 'R_GC%',
            'F_Hairpin_dG_kcal', 'R_Hairpin_dG_kcal',
            'F_Hairpin_Tm', 'R_Hairpin_Tm',
            'F_Homodimer_dG_kcal', 'R_Homodimer_dG_kcal',
            'F_Homodimer_Tm', 'R_Homodimer_Tm',
            'Heterodimer_dG_kcal', 'Heterodimer_Tm',
            'F_BLAST_hits', 'R_BLAST_hits',
            'Nearby_amplicon_count', 'Min_amplicon_size',
            'amplicon_info',
            'QC_Tm_range',
            'QC_Tm_diff',
            'QC_Hairpin',
            'QC_Homodimer',
            'QC_Heterodimer',
            'QC_BLAST_hit',
            'QC_BLAST_amplicon',
            'Final_Result'
        ]) + '\n'
    )

    for line in handle:
        if line.startswith('Forward_Primer'):
            continue

        data = line.strip('\n').split('\t')
        f_name = data[0]
        f_seq = data[1].strip().upper()
        r_name = data[2]
        r_seq = data[3].strip().upper()

        # '+' 포함된 primer는 스킵
        if '+' in f_seq or '+' in r_seq:
            continue

        # ---------- primer3 Thermo 계산 ----------
        # Tm
        f_tm = primer3.calc_tm(f_seq)
        r_tm = primer3.calc_tm(r_seq)

        # GC
        f_gc = gc_fraction(f_seq) * 100.0
        r_gc = gc_fraction(r_seq) * 100.0

        # Hairpin
        f_hp = primer3.calc_hairpin(f_seq)
        r_hp = primer3.calc_hairpin(r_seq)
        f_hp_dg = (f_hp.dg if f_hp.structure_found else 0.0) / 1000.0
        r_hp_dg = (r_hp.dg if r_hp.structure_found else 0.0) / 1000.0
        f_hp_tm = f_hp.tm if f_hp.structure_found else 0.0
        r_hp_tm = r_hp.tm if r_hp.structure_found else 0.0

        # Homodimer
        f_hd = primer3.calc_homodimer(f_seq)
        r_hd = primer3.calc_homodimer(r_seq)
        f_hd_dg = (f_hd.dg if f_hd.structure_found else 0.0) / 1000.0
        r_hd_dg = (r_hd.dg if r_hd.structure_found else 0.0) / 1000.0
        f_hd_tm = f_hd.tm if f_hd.structure_found else 0.0
        r_hd_tm = r_hd.tm if r_hd.structure_found else 0.0

        # Heterodimer
        hetero = primer3.calc_heterodimer(f_seq, r_seq)
        het_dg = (hetero.dg if hetero.structure_found else 0.0) / 1000.0
        het_tm = hetero.tm if hetero.structure_found else 0.0

        # ---------- BLAST ----------

        blast_error = False
        f_hits = r_hits = -1
        nearby_count = -1
        min_amp_size = None
        amp_details = []

        try:
            blast_hits = run_blast_for_primers(f_name, f_seq, r_name, r_seq, blast_db)
            f_hits_list = blast_hits.get(f_name, [])
            r_hits_list = blast_hits.get(r_name, [])
            f_hits = len(f_hits_list)
            r_hits = len(r_hits_list)

            nearby_count, min_amp_size, amp_details = find_nearby_amplicons(
                f_hits_list, r_hits_list, MIN_AMP_BP, MAX_AMP_BP
            )
        except Exception:
            blast_error = True

        # ---------- QC 플래그 (0 / x) ----------

        # 1) Tm 범위
        qc_tm_range = '0' if (TM_MIN <= f_tm <= TM_MAX and TM_MIN <= r_tm <= TM_MAX) else 'x'

        # 2) Tm 차이
        tm_diff = abs(f_tm - r_tm)
        qc_tm_diff = '0' if tm_diff <= TM_DIFF_MAX else 'x'

        # 3) Hairpin
        qc_hairpin = '0' if (f_hp_dg >= HAIRPIN_DG_CUTOFF and r_hp_dg >= HAIRPIN_DG_CUTOFF) else 'x'

        # 4) Homodimer
        qc_homodimer = '0' if (f_hd_dg >= HOMODIMER_DG_CUTOFF and r_hd_dg >= HOMODIMER_DG_CUTOFF) else 'x'

        # 5) Heterodimer
        qc_heterodimer = '0' if het_dg >= HETERODIMER_DG_CUTOFF else 'x'

        # 6) BLAST_hit
        if blast_error:
            qc_blast_hit = 'x'
        else:
            qc_blast_hit = '0' if (f_hits <= BLAST_HIT_MAX and r_hits <= BLAST_HIT_MAX) else 'x'

        # 7) BLAST_amplicon
        if blast_error:
            qc_blast_amplicon = 'x'
        else:
            qc_blast_amplicon = '0' if nearby_count == 0 else 'x'

        # amplicon_info 문자열
        amplicon_info = ';'.join(amp_details) if amp_details else ''

        # 최종 결과
        qc_flags = [
            qc_tm_range, qc_tm_diff, qc_hairpin, qc_homodimer,
            qc_heterodimer, qc_blast_hit, qc_blast_amplicon
        ]
        final_result = 'FAIL' if 'x' in qc_flags else 'PASS'

        # ---------- 출력 ----------
        print(
            '\t'.join([
                f_name, f_seq,
                r_name, r_seq,
                f'{f_tm:.2f}', f'{r_tm:.2f}',
                f'{f_gc:.2f}', f'{r_gc:.2f}',
                f'{f_hp_dg:.2f}', f'{r_hp_dg:.2f}',
                f'{f_hp_tm:.2f}', f'{r_hp_tm:.2f}',
                f'{f_hd_dg:.2f}', f'{r_hd_dg:.2f}',
                f'{f_hd_tm:.2f}', f'{r_hd_tm:.2f}',
                f'{het_dg:.2f}', f'{het_tm:.2f}',
                str(f_hits), str(r_hits),
                str(nearby_count),
                '' if min_amp_size is None else str(min_amp_size),
                amplicon_info,
                qc_tm_range,
                qc_tm_diff,
                qc_hairpin,
                qc_homodimer,
                qc_heterodimer,
                qc_blast_hit,
                qc_blast_amplicon,
                final_result
            ]))
        out.write(
            '\t'.join([
                f_name, f_seq,
                r_name, r_seq,
                f'{f_tm:.2f}', f'{r_tm:.2f}',
                f'{f_gc:.2f}', f'{r_gc:.2f}',
                f'{f_hp_dg:.2f}', f'{r_hp_dg:.2f}',
                f'{f_hp_tm:.2f}', f'{r_hp_tm:.2f}',
                f'{f_hd_dg:.2f}', f'{r_hd_dg:.2f}',
                f'{f_hd_tm:.2f}', f'{r_hd_tm:.2f}',
                f'{het_dg:.2f}', f'{het_tm:.2f}',
                str(f_hits), str(r_hits),
                str(nearby_count),
                '' if min_amp_size is None else str(min_amp_size),
                amplicon_info,
                qc_tm_range,
                qc_tm_diff,
                qc_hairpin,
                qc_homodimer,
                qc_heterodimer,
                qc_blast_hit,
                qc_blast_amplicon,
                final_result
            ]) + '\n'
        )

print(f'완료: 결과 파일 → {output_file}')