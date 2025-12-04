import os
import subprocess

import primer3
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt

# --------- 경로 / 설정 ---------
input_path = '/storage/home/jhkim/Projects/Task/GCX-HearlingLoss_PrimerValidation-2025-12-02/Resources/primer_thermo_result.remove_gc_clamp.tsv'
output_dir = '/storage/home/jhkim/Projects/Task/GCX-HearlingLoss_PrimerValidation-2025-12-02/Results'
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'primer_thermo_blast_result.tab')

# BLAST 실행 파일 / DB
BLASTN = "/storage/home/jhkim/Apps/ncbi-blast-2.17.0+/bin/blastn"
BLASTDBCMD = "/storage/home/jhkim/Apps/ncbi-blast-2.17.0+/bin/blastdbcmd"
blast_db = '/storage/home/jhkim/Apps/ncbi-blast-2.17.0+/ref/hg38'

# BLAST 필터 기준
BLAST_IDENTITY_THRESHOLD = 80.0   # %
BLAST_LENGTH_THRESHOLD = 10       # 최소 매칭 길이 (nt)
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
BLAST_HIT_MAX = 1  # on-target 1개만 허용

# F/R 같은 chr 내 잠재 amplicon 거리 기준 (bp)
MIN_AMP_BP = 50
MAX_AMP_BP = 300

# self-amplicon 그림 저장 폴더
self_plot_dir = os.path.join(output_dir, "self_ascii_plots")
os.makedirs(self_plot_dir, exist_ok=True)


# ---------- 유틸 함수들 ----------

def get_reference_subseq(db, chrom, start, end):
    """
    blastdbcmd로 ref DB에서 특정 위치 서열 가져오기.
    start, end: 1-based inclusive 좌표
    """
    cmd = [
        BLASTDBCMD,
        "-db", db,
        "-entry", chrom,
        "-range", f"{start}-{end}",
        "-strand", "plus",
        "-outfmt", "%s"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    result.check_returncode()
    seq = "".join(result.stdout.splitlines()).upper()
    return seq


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


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


def run_blast_for_primers(f_name, f_seq, r_name, r_seq, db):
    """
    Forward, Reverse 프라이머를 하나의 FASTA로 만들어 blastn 수행 후
    각 primer별 hit 목록을 반환.
    hit dict에 qseq/sseq까지 포함.
    """
    orig_f_name = f_name  # BLAST qseqid용
    orig_r_name = r_name

    # 파일 이름용 safe name
    safe_f_name = (
        orig_f_name.replace('(', '_')
        .replace(')', '_')
        .replace(' ', '_')
        .replace('/', '_')
    )
    fasta_str = f">{orig_f_name}\n{f_seq}\n>{orig_r_name}\n{r_seq}\n"

    fasta_path = f'{safe_f_name}.fasta'

    with open(fasta_path, 'w') as tmp_fasta:
        tmp_fasta.write(fasta_str)

    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"

    cmd = [
        BLASTN,
        "-task", "blastn-short",
        "-db", db,
        "-query", fasta_path,
        "-outfmt", outfmt,
        "-num_alignments", str(BLAST_MAX_ALIGNMENTS),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        # BLAST 에러 시 빈 hit 반환
        return {orig_f_name: [], orig_r_name: []}

    hits = {
        orig_f_name: [],
        orig_r_name: [],
    }

    for line in result.stdout.strip().splitlines():
        if not line.strip():
            continue
        cols = line.strip('\n').split("\t")
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

        # 필터 기준 적용
        if pident < BLAST_IDENTITY_THRESHOLD or length < BLAST_LENGTH_THRESHOLD:
            continue

        if qseqid not in hits:
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
            "qseq": qseq,
            "sseq": sseq,
        })

    return hits


def make_primer_full_annotation(hit, primer_seq):
    """
    primer 전체 길이 기준으로:
      - primer_full: mismatch 위치는 소문자
      - match_mask : '|' = match, 'x' = mismatch, ' ' = 미정렬
    를 만들어 반환한다.
    BLAST qstart/qseq/sseq 기준 (갭 없다는 가정, short oligo에서 보통 그럼)
    """
    qseq_aln = hit["qseq"]
    sseq_aln = hit["sseq"]
    qstart = hit["qstart"]  # 1-based

    plen = len(primer_seq)
    primer_chars = list(primer_seq)   # 기본은 그대로
    mask_chars = [' '] * plen

    qpos = qstart - 1  # 0-based index in primer_seq

    for qb, sb in zip(qseq_aln, sseq_aln):
        if qb == '-':
            # query gap → primer index 안 움직임
            continue
        if qpos < 0 or qpos >= plen:
            break

        match = (qb.upper() == sb.upper())
        if match:
            mask_chars[qpos] = '|'
        else:
            mask_chars[qpos] = 'x'
            # mismatch 위치 primer는 소문자로
            primer_chars[qpos] = primer_chars[qpos].lower()
        qpos += 1

    primer_full = ''.join(primer_chars)
    mask_line = ''.join(mask_chars)
    return primer_full, mask_line


def plot_self_amplicon_ascii_style(
    db,
    primer_name,
    primer_seq,
    chrom,
    hit_i,
    hit_j,
    label,      # 'F-SELF' or 'R-SELF'
    out_png
):
    """
    ref 한 줄 + primer1/primer2를 ref 좌표에 맞게 정렬해서
    ASCII 스타일로 그림(PNG) 저장.
    각 primer에 대해:
      - 전체 primer 서열 (mismatch는 소문자)
      - match_mask ('|' / 'x' / ' ')
    도 함께 그려 줌.
    """

    def span_coords(h):
        a = h["sstart"]
        b = h["send"]
        return min(a, b), max(a, b)

    span_i_start, span_i_end = span_coords(hit_i)
    span_j_start, span_j_end = span_coords(hit_j)

    region_start = min(span_i_start, span_j_start)
    region_end   = max(span_i_end, span_j_end)

    ref_seq = get_reference_subseq(db, chrom, region_start, region_end)
    L = len(ref_seq)

    # --- ref 위에 찍을 primer 라인 (match 대문자 / mismatch 소문자) ---
    def build_ref_aligned_line(hit, region_start, primer_seq):
        strand, _3p = hit_strand_and_3end(hit)
        sstart = hit["sstart"]
        send   = hit["send"]
        qseq_aln = hit["qseq"]
        sseq_aln = hit["sseq"]

        # ref 기준 시작 위치
        # (subject alignment 상에서 sseq_aln의 첫 non-gap이 region_start 기준 어디인지 추적)
        span_start = min(sstart, send)
        offset = span_start - region_start  # 대략적인 시작 offset

        line_chars = [' '] * L

        # 간단 버전: 갭 거의 없다고 가정하고, offset부터 한 글자씩 채움
        # match → primer base 대문자, mismatch → primer base 소문자
        idx_ref = offset
        for qb, sb in zip(qseq_aln, sseq_aln):
            if qb == '-':
                continue
            if idx_ref < 0 or idx_ref >= L:
                break
            match = (qb.upper() == sb.upper())
            base = qb.upper() if match else qb.lower()
            line_chars[idx_ref] = base
            idx_ref += 1

        return ''.join(line_chars)

    primer1_ref_line = build_ref_aligned_line(hit_i, region_start, primer_seq)
    primer2_ref_line = build_ref_aligned_line(hit_j, region_start, primer_seq)

    # --- primer 전체 서열 기준 match/mismatch annotation ---
    primer_full_1, mask_1 = make_primer_full_annotation(hit_i, primer_seq)
    primer_full_2, mask_2 = make_primer_full_annotation(hit_j, primer_seq)

    # --- 그림 내용 구성 ---
    title = f"{label} for {primer_name} on {chrom}:{region_start}-{region_end}"

    # 그림 그리기
    fig, ax = plt.subplots(figsize=(min(18, L / 4 + 6), 5))
    ax.axis('off')

    y = 0.94
    ax.text(0.01, y, title, fontsize=11, fontfamily="monospace", transform=ax.transAxes)
    y -= 0.08

    # ref + 두 primer (ref 좌표에 맞춰 정렬)
    ax.text(0.01, y, f"{'ref':8s} {ref_seq}", fontsize=10, fontfamily="monospace", transform=ax.transAxes)
    y -= 0.08
    ax.text(0.01, y, f"{'primer1':8s} {primer1_ref_line}", fontsize=10, fontfamily="monospace", transform=ax.transAxes)
    y -= 0.08
    ax.text(0.01, y, f"{'primer2':8s} {primer2_ref_line}", fontsize=10, fontfamily="monospace", transform=ax.transAxes)
    y -= 0.10

    # primer 전체 서열 + match/mismatch mask (site1)
    ax.text(0.01, y, "Site1 (primer full vs alignment):", fontsize=10,
            fontfamily="monospace", transform=ax.transAxes)
    y -= 0.06
    ax.text(0.01, y, f"  primer1_full: {primer_full_1}", fontsize=10,
            fontfamily="monospace", transform=ax.transAxes)
    y -= 0.06
    ax.text(0.01, y, f"  match_mask  : {mask_1}", fontsize=10,
            fontfamily="monospace", transform=ax.transAxes)
    y -= 0.10

    # primer 전체 서열 + match/mismatch mask (site2)
    ax.text(0.01, y, "Site2 (primer full vs alignment):", fontsize=10,
            fontfamily="monospace", transform=ax.transAxes)
    y -= 0.06
    ax.text(0.01, y, f"  primer2_full: {primer_full_2}", fontsize=10,
            fontfamily="monospace", transform=ax.transAxes)
    y -= 0.06
    ax.text(0.01, y, f"  match_mask  : {mask_2}", fontsize=10,
            fontfamily="monospace", transform=ax.transAxes)

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def find_nearby_amplicons(f_hits, r_hits,
                          min_bp=MIN_AMP_BP,
                          max_bp=MAX_AMP_BP,
                          f_len=None,
                          r_len=None):
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
            if f_strand == '+' and r_strand == '-' and f_3p < r_3p:
                valid_orientation = True
            elif f_strand == '-' and r_strand == '+' and r_3p < f_3p:
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


def find_self_amplicons(hits,
                        min_bp=MIN_AMP_BP,
                        max_bp=MAX_AMP_BP,
                        primer_len=None,
                        label='F',
                        primer_name='PRIMER',
                        primer_seq=None,
                        db=None,
                        plot_dir=None):
    """
    같은 primer(hits) 안에서 self-amplicon 찾기 + (옵션) ref+primer 그림 저장

    - 같은 chr
    - 서로 반대 strand
    - 서로를 향하는 방향
    - 길이: core distance + 2 * primer_len

    label: 'F' 또는 'R' (F-SELF / R-SELF 구분용)
    primer_seq: 전체 primer 서열(5'→3')
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
            if strand_i == '+' and strand_j == '-' and p3_i < p3_j:
                valid_orientation = True
            elif strand_i == '-' and strand_j == '+' and p3_j < p3_i:
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

                # 그림 저장 (옵션)
                if db is not None and plot_dir is not None:
                    png_name = f"{primer_name}_{label}-SELF_{chr_i}_{p3_i}_{p3_j}.png"
                    out_png = os.path.join(plot_dir, png_name)
                    plot_self_amplicon_ascii_style(
                        db=db,
                        primer_name=primer_name,
                        primer_seq=primer_seq,
                        chrom=chr_i,
                        hit_i=hi,
                        hit_j=hj,
                        label=f"{label}-SELF",
                        out_png=out_png
                    )

    return count, min_size, details


# ---------- 메인 루프 ----------

with open(input_path, 'r') as handle, open(output_file, 'w') as out:

    # 결과 헤더
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
        if len(data) < 4:
            continue

        f_name = data[0]
        f_seq = data[1].strip().upper()
        r_name = data[2]
        r_seq = data[3].strip().upper()

        # '+' 포함된 primer는 스킵
        if '+' in f_seq or '+' in r_seq:
            continue

        # ---------- primer3 Thermo 계산 ----------
        f_tm = primer3.calc_tm(f_seq)
        r_tm = primer3.calc_tm(r_seq)

        f_gc = gc_fraction(f_seq) * 100.0
        r_gc = gc_fraction(r_seq) * 100.0

        f_hp = primer3.calc_hairpin(f_seq)
        r_hp = primer3.calc_hairpin(r_seq)
        f_hp_dg = (f_hp.dg if f_hp.structure_found else 0.0) / 1000.0
        r_hp_dg = (r_hp.dg if r_hp.structure_found else 0.0) / 1000.0
        f_hp_tm = f_hp.tm if f_hp.structure_found else 0.0
        r_hp_tm = r_hp.tm if r_hp.structure_found else 0.0

        f_hd = primer3.calc_homodimer(f_seq)
        r_hd = primer3.calc_homodimer(r_seq)
        f_hd_dg = (f_hd.dg if f_hd.structure_found else 0.0) / 1000.0
        r_hd_dg = (r_hd.dg if r_hd.structure_found else 0.0) / 1000.0
        f_hd_tm = f_hd.tm if f_hd.structure_found else 0.0
        r_hd_tm = r_hd.tm if r_hd.structure_found else 0.0

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

            # F-R amplicon
            c_FR, min_FR, det_FR = find_nearby_amplicons(
                f_hits_list, r_hits_list,
                MIN_AMP_BP, MAX_AMP_BP,
                f_len=len(f_seq), r_len=len(r_seq)
            )
            # F-SELF
            c_FF, min_FF, det_FF = find_self_amplicons(
                f_hits_list,
                MIN_AMP_BP, MAX_AMP_BP,
                primer_len=len(f_seq),
                label='F',
                primer_name=f_name,
                primer_seq=f_seq,
                db=blast_db,
                plot_dir=self_plot_dir
            )
            # R-SELF
            c_RR, min_RR, det_RR = find_self_amplicons(
                r_hits_list,
                MIN_AMP_BP, MAX_AMP_BP,
                primer_len=len(r_seq),
                label='R',
                primer_name=r_name,
                primer_seq=r_seq,
                db=blast_db,
                plot_dir=self_plot_dir
            )

            nearby_count = c_FR + c_FF + c_RR
            mins = [x for x in [min_FR, min_FF, min_RR] if x is not None]
            min_amp_size = min(mins) if mins else None
            amp_details = det_FR + det_FF + det_RR

        except Exception:
            blast_error = True

        # ---------- QC (O / X) ----------
        qc_tm_range = 'O' if (TM_MIN <= f_tm <= TM_MAX and TM_MIN <= r_tm <= TM_MAX) else 'X'

        tm_diff = abs(f_tm - r_tm)
        qc_tm_diff = 'O' if tm_diff <= TM_DIFF_MAX else 'X'

        qc_hairpin = 'O' if (f_hp_dg >= HAIRPIN_DG_CUTOFF and r_hp_dg >= HAIRPIN_DG_CUTOFF) else 'X'

        qc_homodimer = 'O' if (f_hd_dg >= HOMODIMER_DG_CUTOFF and r_hd_dg >= HOMODIMER_DG_CUTOFF) else 'X'

        qc_heterodimer = 'O' if het_dg >= HETERODIMER_DG_CUTOFF else 'X'

        if blast_error:
            qc_blast_hit = 'X'
            qc_blast_amplicon = 'X'
        else:
            qc_blast_hit = 'O' if (f_hits <= BLAST_HIT_MAX and r_hits <= BLAST_HIT_MAX) else 'X'
            # amplicon 하나라도 있으면 FAIL
            qc_blast_amplicon = 'O' if nearby_count == 0 else 'X'

        amplicon_info = ';'.join(amp_details) if amp_details else ''

        # 최종 결과는 BLAST_amplicon 포함 QC 기준으로
        qc_flags = [
            qc_tm_range, qc_tm_diff, qc_hairpin, qc_homodimer,
            qc_heterodimer, qc_blast_amplicon
        ]
        final_result = 'FAIL' if 'X' in qc_flags else 'PASS'

        row = [
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
        ]

        print('\t'.join(row))
        out.write('\t'.join(row) + '\n')

print(f'완료: 결과 파일 → {output_file}')
print(f'SELF amplicon 그림 폴더 → {self_plot_dir}')
