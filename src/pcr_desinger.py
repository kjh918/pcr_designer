import pysam
import pandas as pd
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Optional
from Bio.Seq import Seq

def revcomp(seq: str) -> str:
    """Reverse-complement a DNA sequence."""
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]

def build_pair_alignment_lines(
        ref_full: str,
        primer_seq: str,
        offset: int,
        primer_label: str,
        ref_label: str,
    ):
    """
    ref_full: 그룹 전체 구간의 reference 서열 (한 줄)
    primer_seq: 5'->3' primer 서열 (ref_full과 같은 방향으로 정렬된 것)
    offset: primer가 ref_full 어디서 시작하는지 (0-based index)
    primer_label, ref_label: 출력시 앞에 붙일 라벨
    """
    L = len(ref_full)

    # primer / match 라인용 배열 (처음엔 공백)
    arr_primer = [" "] * L
    arr_match = [" "] * L

    for i, base in enumerate(primer_seq):
        pos = offset + i
        if 0 <= pos < L:
            arr_primer[pos] = base
            if base == ref_full[pos]:
                arr_match[pos] = "|"

    line_primer = f"{primer_label}: " + "".join(arr_primer)
    line_match  = " " * (len(primer_label) + 2) + "".join(arr_match)
    line_ref    = f"{ref_label}: " + ref_full

    return line_primer, line_match, line_ref

@dataclass
class PrimerRecord:
    chrom: str
    start: int
    end: int
    name: str
    sequence: str      # 입력 primer 서열 (5'->3')
    strand: str        # '+', '-'
    ref_seq: Optional[str] = None          # reference에서 뽑은 서열 (ref 방향)
    bind_seq: Optional[str] = None         # reference에 실제 붙는 primer 서열(ref 방향)
    matches: Optional[bool] = None         # 전체 일치 여부
    per_base_match: Optional[List[bool]] = None  # 각 위치별 match / mismatch

class PrimerBindingAligner:
    """
    BED + reference FASTA 기준으로
    primer binding 서열을 group 단위로 관리하고
    정렬 형태로 시각화하는 클래스.
    """

    def __init__(self, fasta_path: str, bed_path: str):
        self.fasta_path = fasta_path
        self.bed_path = bed_path
        self.ref = pysam.FastaFile(fasta_path)

        self.primers: List[PrimerRecord] = []
        self.groups: dict[str, List[PrimerRecord]] = {}  # 추가됨!
        self.ref_templete: dict[str, List[PrimerRecord]] = {} 
        self.group_region: dict[str, dict] = {}     

    def load_bed(
        self,
        chrom_col: int = 0,
        start_col: int = 1,
        end_col: int = 2,
        name_col: int = 3,
        seq_col: int = 4,
        strand_col: int = 5,
        group_col: int | None = None,   # ★ 그룹 컬럼 추가 (선택)
        header: bool = True,
        sep: str = "\t",
    ):
        """
        BED/TSV에서 primer 읽고 group별로 dict에 저장.
        
        기본 포맷: chrom start end name sequence strand (group 컬럼은 선택)
        """
        df = pd.read_csv(
            self.bed_path,
            sep=sep,
            header=None,
            dtype={chrom_col: str},
        )


        for _, row in df.iterrows():
            chrom = row.iloc[0]
            start = int(row.iloc[1])
            end = int(row.iloc[2])
            name = str(row.iloc[3])
            seq = str(row.iloc[4]).replace(" ", "").upper()
            strand = str(row.iloc[5])
            group = str(row.iloc[6])

            primer = PrimerRecord(
                chrom=chrom,
                start=start,
                end=end,
                name=name,
                sequence=seq,
                strand=strand,
            )

            # 전체 리스트에도 저장
            self.primers.append(primer)

            # 그룹 dict에도 저장
            if group not in self.groups:
                self.groups[group] = []
            self.groups[group].append(primer)

    def check():
        pass 
    def compute_binding(self, extend=30, reverse=False):
        """
        group별 reference 범위를 확정하고,
        - ref_seq (lower)
        - target / mismatch 만 UPPER
        - reverse_complement ref_seq도 생성
        """

        for group_id, primer_list in self.groups.items():

            chrom = ''
            target_start = None
            target_end   = None

            # target / mismatch 위치 저장
            highlight_positions = []  # 리스트에 (start,end,type) 저장

            for p in primer_list:
                if p.name == 'primer_1':
                    chrom = p.chrom
                    target_start = p.start - extend 
                    primer1 = p.sequence

                elif p.name == 'primer_2':
                    target_end = p.end + extend
                    primer2 = p.sequence[::-1]

                # target / mismatch 영역 수집
                if p.name == "target" or p.name == "missmatch":
                    highlight_positions.append((p.start, p.end))

            # fallback (플로우 보호)
            if chrom == '':
                chrom = primer_list[0].chrom
            if target_start is None:
                target_start = min(p.start for p in primer_list) - extend
            if target_end is None:
                target_end = max(p.end for p in primer_list) + extend
        
            # reference 가져오기 (lowercase)
            ref_seq = self.ref.fetch(chrom, target_start, target_end).lower()

            # --------------------------
            # ⭐ target / mismatch 영역만 UPPER로 변환
            # --------------------------
            # ref_list = list(ref_seq)

            # for (hs, he) in highlight_positions:
            #     local_s = hs - target_start
            #     local_e = he - target_start
                
            #     # ref_seq 길이 안쪽으로만 자르기
            #     local_s = max(0, local_s)
            #     local_e = min(len(ref_list), local_e)

            #     for i in range(local_s, local_e):
            #         ref_list[i] = ref_list[i].upper()

            # ref_seq = "".join(ref_list)

            # reverse complement (그대로)
            reverse_ref_seq = str(Seq(ref_seq).complement())

            # 저장
            self.ref_templete[group_id] = [ref_seq, reverse_ref_seq]
            # primer_1 시작 위치 (0-based index on ref_seq)
            # print 확인
            print(f"\n=== GROUP {group_id} ===")
            print(f"=== REF       hg38 ===")
            print("forward  :", ' '*(extend) + primer1)
            print("          ", ' '*(extend) + '|'*len(primer1))
            print("REF(5->3):", ref_seq.upper())
            print("          ", '+'*len(ref_seq))
            print("REF(3->5):", reverse_ref_seq.upper())
            print("          ", ' '*(primer_list[3].start - target_start) + '|'*len(primer2))
            print("reverse  :", ' '*(primer_list[3].start - target_start) + primer2)






# 사용 예시
if __name__ == "__main__":


    fasta = "/storage/home/jhkim/Projects/Task/GCX-PCRPrimerDesign-2025-11-27/Resources/Reference/hg38/hg38.fa"   # reference FASTA
    bed = "/storage/home/jhkim/Projects/Task/GCX-PCRPrimerDesign-2025-11-27/Resources/Bed/target_primer.test.bed"      # primer BED/TSV

    aligner = PrimerBindingAligner(fasta, bed)
    aligner.load_bed()        # chrom start end name sequence strand 기본 포맷 가정
    aligner.compute_binding()


    # group별 ref/primer 정렬 출력
    # for gid in aligner.groups.keys():
    #     aligner.print_group_pair_alignments(gid)

    # 3) 특정 primer 하나 mismatch 빨간색으로 그림
    # aligner.plot_alignment("primer1")
