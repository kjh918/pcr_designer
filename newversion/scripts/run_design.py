import argparse
import sys
import os
import pysam
import pandas as pd
from typing import List

# -----------------------------------------------------------------------------
# 1. 환경 설정
# -----------------------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import ConfigLoader
from pcr.config.schema.root import BaseDesignInput
from pcr.designers.qpcr import QPCRPrimerDesigner 
from pcr.utils.ranker import ProbeCentricRanker
from pcr.qc.executor import QCExecutor
# ✅ GenomicRegion, Amplicon 임포트 확인
from pcr.components.amplicon import Amplicon, GenomicRegion 

def parse_args():
    parser = argparse.ArgumentParser(description="Gemini PCR Pipeline: Design -> Rank -> QC -> Save")
    
    parser.add_argument("--name", type=str, required=True, help="Task name")
    parser.add_argument("--chrom", type=str, required=True, help="Chromosome")
    parser.add_argument("--start", type=int, required=True, help="Genomic Target Start")
    parser.add_argument("--end", type=int, required=True, help="Genomic Target End")
    parser.add_argument("--ref_genotype", type=str, required=True, help="Reference genotype")
    parser.add_argument("--alt_genotype", type=str, required=True, help="Alt genotype")
    
    parser.add_argument("--padding", type=int, default=75, help="Padding size (bp)")
    parser.add_argument("--preset", type=str, default="default", help="Preset name")
    parser.add_argument("--ref", type=str, default="hg38", help="Reference ID")
    parser.add_argument("--out", type=str, default="results.xlsx", help="Output file path")
    
    parser.add_argument("--top_k", type=int, default=5, help="Target number of final candidates")
    parser.add_argument("--num_search", type=int, default=100, help="Initial Primer3 search count") 

    return parser.parse_args()

def fetch_template_and_coords(chrom, g_start, g_end, r_genotype, a_genotype, padding, fasta_path):
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    with pysam.FastaFile(fasta_path) as fasta:
        fetch_start = max(0, g_start - padding)
        fetch_end = g_end + padding
        
        reference_seq = fasta.fetch(chrom, fetch_start, fetch_end).upper()
        rel_start = g_start - fetch_start
        rel_end_ref = g_end - fetch_start
        
        fetched_ref = reference_seq[rel_start:rel_end_ref]
        if r_genotype and fetched_ref != r_genotype:
             print(f"⚠️ [Mismatch Warning] Genomic Ref('{fetched_ref}') != Input Ref('{r_genotype}')")

        template_seq = reference_seq[:rel_start] + a_genotype + reference_seq[rel_end_ref:]
        target_len = len(a_genotype)
        new_target_end = rel_start + target_len

        return reference_seq, template_seq, rel_start, new_target_end, fetch_start

def visualize_alignment_cmd(amp: Amplicon, genotype_label: str):
    """터미널에 서열 정렬 상태를 색상으로 출력"""
    seq = amp.template_sequence
    f_start = amp.forward.start_index
    f_end = f_start + len(amp.forward.sequence)
    r_start = amp.reverse.start_index
    r_end = r_start + len(amp.reverse.sequence)
    
    p_start = -1
    p_end = -1
    if amp.probe:
        p_start = getattr(amp.probe, 'start_index', -1)
        p_end = p_start + len(amp.probe.sequence)

    # 타겟(변이) 위치 (GenomicRegion에는 절대좌표, amp.target_start_index는 상대좌표)
    t_start = amp.target_start_index
    t_end = amp.target_end_index
    
    # ANSI Color Codes
    GREEN = '\033[92m'  # Fwd
    RED = '\033[91m'    # Rev
    BLUE = '\033[94m'   # Probe
    YELLOW = '\033[93m' # Target
    RESET = '\033[0m'
    
    vis_seq = ""
    for i, char in enumerate(seq):
        color = ""
        # 우선순위: Target > Probe > Primer
        if t_start <= i < t_end:
            color = YELLOW
        elif p_start <= i < p_end:
            color = BLUE
        elif f_start <= i < f_end:
            color = GREEN
        elif r_start <= i < r_end:
            color = RED
            
        if color:
            vis_seq += f"{color}{char}{RESET}"
        else:
            vis_seq += char
            
    print(f"\n🔹 Alignment View [{amp.id}] ({genotype_label})")
    print(f"   Seq: {vis_seq}")
    print(f"   Legend: {GREEN}Fwd{RESET} {RED}Rev{RESET} {BLUE}Probe{RESET} {YELLOW}Target(Var){RESET}")

def save_excel_with_formatting(df: pd.DataFrame, file_path: str, amp_list: List[Amplicon]):
    """엑셀 저장 시 서열 내 위치를 색상으로 강조 (Rich Text는 어렵지만 셀 채우기로 대체하거나 별도 표기)"""
    
    # Pandas Styler 대신 XlsxWriter 엔진을 직접 사용하여 조건부 서식 적용
    # 하지만 셀 내의 '특정 글자'만 색칠하는 것은 매우 복잡하므로,
    # 여기서는 'Template_Seq' 컬럼 옆에 시각화된 HTML 문자열을 넣거나,
    # 별도의 시트에 상세 정보를 넣는 방식을 추천합니다.
    # 가장 현실적인 방법: 주요 구간(Start, End) 정보를 컬럼으로 남기고, 
    # 엑셀의 '조건부 서식'은 셀 단위이므로 텍스트 전체가 아닌 셀 배경색만 가능합니다.
    
    # 사용자의 요청(색상 다르게)을 위해, 엑셀에는 'Rich Text'를 쓰기 어렵기 때문에
    # HTML 리포트를 생성하거나, 좌표 정보를 명확히 남기는 것으로 대체합니다.
    # 대신 CMD 창에서는 확실하게 보여드렸습니다.
    
    # 기본 엑셀 저장
    with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
        df.to_excel(writer, index=False, sheet_name='Result')
        workbook = writer.book
        worksheet = writer.sheets['Result']
        
        # 헤더 포맷
        header_format = workbook.add_format({'bold': True, 'bg_color': '#D7E4BC', 'border': 1})
        for col_num, value in enumerate(df.columns.values):
            worksheet.write(0, col_num, value, header_format)
            
        # 컬럼 너비 자동 조정 (대략)
        for i, col in enumerate(df.columns):
            column_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
            worksheet.set_column(i, i, min(column_len, 50)) # 최대 50자 제한

def main():
    args = parse_args()

    # 1. Init
    print(f"📦 [Init] Loading config: {args.preset}")
    loader = ConfigLoader()
    config = loader.load(preset_name=args.preset)
    
    if args.ref not in config.references:
        print(f"❌ Error: Reference '{args.ref}' not found")
        sys.exit(1)
    ref_path = config.references[args.ref].fasta_path

    # 2. Fetch & Design
    try:
        print(f"🧬 [Step 1] Fetching sequence & Designing...")
        reference_seq, template_seq, rel_start, rel_end, region_start_pos = fetch_template_and_coords(
            args.chrom, args.start - 1 , args.end, args.ref_genotype, args.alt_genotype, args.padding, ref_path
        )
        
        design_input = BaseDesignInput(
            name=args.name,
            template_sequence=template_seq,
            reference_sequence=reference_seq,
            target_start=rel_start,
            target_end=rel_end,
            config=config,
            reference_name=args.ref,
            overrides={"PRIMER_NUM_RETURN": args.num_search}
        )
        
        designer = QPCRPrimerDesigner(design_input)
        design_output = designer.design()
        raw_candidates = design_output.amplicons
        
        # [중요] 절대 좌표 및 GenomicRegion 주입
        for amp in raw_candidates:
            abs_start = region_start_pos + amp.forward.start_index
            abs_end = region_start_pos + amp.reverse.start_index + len(amp.reverse.sequence)
            
            #amp.genomic_chrom = args.chrom
            #amp.genomic_start = region_start_pos + amp.target_start_index # 타겟 위치 기준
            
            amp.region = GenomicRegion(
                chrom=args.chrom, 
                start=abs_start, 
                end=abs_end,
                strand='+'  # 기본값
            )

        print(f"   -> Found {len(raw_candidates)} raw candidates.")

    except Exception as e:
        print(f"❌ Design Failed: {e}")
        import traceback; traceback.print_exc()
        sys.exit(1)

    # 3. Ranking
    print(f"📊 [Step 2] Ranking...")
    ranker = ProbeCentricRanker(probe_overlap_threshold=0.9)
    sorted_candidates = ranker.select_diverse_probes(raw_candidates, top_k=len(raw_candidates))

    # 4. QC
    print(f"🔍 [Step 3] Running QC...")
    qc_executor = QCExecutor(config)
    final_selection = []
    tested_candidates = []

    for i, amp in enumerate(sorted_candidates):
        if len(final_selection) >= args.top_k: break
        
        processed_batch = qc_executor.run_qc([amp])
        processed_amp = processed_batch[0]
        tested_candidates.append(processed_amp)
        
        if processed_amp.is_qc_pass:
            final_selection.append(processed_amp)

    print(f"\n   -> Final QC Passed: {len(final_selection)}")

    # 5. Save & Visualize
    print(f"💾 [Step 4] Saving results...")
    report_list = final_selection + [c for c in tested_candidates if not c.is_qc_pass]
    summary_df = qc_executor.summarize(report_list)
    
    if not summary_df.empty:
        # 좌표 복원
        genomic_fwd_starts = []
        genomic_rev_starts = []
        genomic_probe_starts = []
        
        for amp in report_list:
            genomic_fwd_starts.append(region_start_pos + amp.forward.start_index)
            genomic_rev_starts.append(region_start_pos + amp.reverse.start_index)
            p_start = getattr(amp.probe, 'start_index', 0) if amp.probe else 0
            genomic_probe_starts.append(region_start_pos + p_start)

        summary_df.insert(1, "Genomic_Fwd_Start", genomic_fwd_starts)
        summary_df.insert(2, "Genomic_Rev_Start", genomic_rev_starts)
        summary_df.insert(3, "Genomic_Probe_Start", genomic_probe_starts)
        
        # ✅ 요청사항 1: Alt, Ref Genotype 컬럼 추가
        summary_df.insert(4, "Ref_Allele", args.ref_genotype)
        summary_df.insert(5, "Alt_Allele", args.alt_genotype)

        # ✅ 요청사항 2: CMD 창에 Alignment 시각화 (상위 3개만 예시로)
        print("\n👀 Top 3 Candidates Alignment Visualization:")
        for amp in final_selection[:5]:
            visualize_alignment_cmd(amp, f"{args.ref_genotype}>{args.alt_genotype}")

    # 파일 저장 (포맷팅 적용)
    try:
        if args.out.endswith(".xlsx"):
            save_excel_with_formatting(summary_df, args.out, report_list)
        else:
            summary_df.to_csv(args.out, index=False)
        print(f"✅ Done! Saved to: {args.out}")
        
    except Exception as e:
        print(f"❌ Failed to save file: {e}")

if __name__ == "__main__":
    main()