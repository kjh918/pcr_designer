import argparse
import sys
import os
import pysam
import pandas as pd

# 경로 설정
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import ConfigLoader
from pcr.config.schema.root import BaseDesignInput
from pcr.designers.base import BasePrimerDesigner
from pcr.utils.ranker import ProbeCentricRanker
from pcr.qc.executor import QCExecutor

def parse_args():
    parser = argparse.ArgumentParser(description="Gemini PCR Pipeline: Design -> QC -> Rank -> Save")
    
    parser.add_argument("--name", type=str, required=True, help="Task name")
    parser.add_argument("--chrom", type=str, required=True, help="Chromosome (e.g. chr7)")
    parser.add_argument("--start", type=int, required=True, help="Genomic Target Start")
    parser.add_argument("--end", type=int, required=True, help="Genomic Target End")
    
    parser.add_argument("--padding", type=int, default=300, help="Padding size (bp)")
    parser.add_argument("--preset", type=str, default="default", help="Preset name")
    parser.add_argument("--ref", type=str, default="hg38", help="Reference ID")
    parser.add_argument("--out", type=str, default="results.xlsx", help="Output file")
    
    parser.add_argument("--top_k", type=int, default=10, help="Final selection count")
    parser.add_argument("--num_search", type=int, default=10, help="Primer3 raw count") 

    return parser.parse_args()

def fetch_template_and_coords(chrom, g_start, g_end, padding, fasta_path):
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")
    
    with pysam.FastaFile(fasta_path) as fasta:
        fetch_start = max(0, g_start - padding)
        fetch_end = g_end + padding
        
        # 서열 추출 (대문자 변환)
        template_seq = fasta.fetch(chrom, fetch_start, fetch_end).upper()
        
        # 상대 좌표 계산
        rel_start = g_start - fetch_start
        rel_end = rel_start + (g_end - g_start)
        
        return template_seq, rel_start, rel_end, fetch_start

def main():
    args = parse_args()

    # -------------------------------------------------------------
    # 0. 설정 및 데이터 로드
    # -------------------------------------------------------------
    print(f"📦 [Init] Loading config: {args.preset}")
    loader = ConfigLoader()
    config = loader.load(preset_name=args.preset)

    if args.ref not in config.references:
        print(f"❌ Error: Reference '{args.ref}' not found in system.yaml")
        sys.exit(1)
    
    ref_path = config.references[args.ref].fasta_path

    try:
        print(f"🧬 [Fetch] Fetching sequence {args.chrom}:{args.start}-{args.end}...")
        template_seq, rel_start, rel_end, region_start_pos = fetch_template_and_coords(
            args.chrom, args.start, args.end, args.padding, ref_path
        )
    except Exception as e:
        print(f"❌ Fetch Failed: {e}")
        sys.exit(1)

    design_input = BaseDesignInput(
        name=args.name,
        template_sequence=template_seq,
        target_start=rel_start,
        target_end=rel_end,
        config=config,
        reference_name=args.ref,
        overrides={"PRIMER_NUM_RETURN": args.num_search}
    )

    # -------------------------------------------------------------
    # 1. 디자인 (Design) - Primer3 실행
    # -------------------------------------------------------------
    print(f"🚀 [Step 1] Running Primer3 (Target: {args.num_search})...")
    designer = BasePrimerDesigner(design_input)
    design_output = designer.design()

    if not design_output.amplicons:
        print("❌ Primer3 found 0 candidates.")
        sys.exit(1)
    
    raw_candidates = design_output.amplicons
    print(f"   -> Found {len(raw_candidates)} raw candidates.")

    # -------------------------------------------------------------
    # 2. 검증 (QC) - 모든 후보에 대해 수행
    # -------------------------------------------------------------
    print(f"🔍 [Step 2] Running QC on ALL candidates...")
    qc_executor = QCExecutor(config)
    qc_processed = qc_executor.run_qc(raw_candidates)

    # QC 통과한 것만 필터링
    passed_candidates = [amp for amp in qc_processed if amp.is_qc_pass]
    print(f"   -> QC Passed: {len(passed_candidates)} / {len(raw_candidates)}")

    if not passed_candidates:
        print("⚠️ Warning: No candidates passed QC. Saving failed results.")
        # 실패했어도 로그를 위해 전체 저장할지 결정
        # 여기서는 실패한 것들을 그대로 둠
        final_selection = qc_processed 
    else:
        # -------------------------------------------------------------
        # 3. 랭킹 (Ranking) - 통과한 것 중에서 Top K 선정
        # -------------------------------------------------------------
        print(f"📊 [Step 3] Ranking passed candidates (Diversity check)...")
        ranker = ProbeCentricRanker(probe_overlap_threshold=0.9)
        
        # QC 통과한 것들 중에서만 랭킹 산정
        final_selection = ranker.select_diverse_probes(
            passed_candidates, 
            top_k=args.top_k
        )
        print(f"   -> Selected Top {len(final_selection)} candidates.")

    # -------------------------------------------------------------
    # 4. 결과 저장 (Save) - 절대 좌표 변환 포함
    # -------------------------------------------------------------
    print(f"💾 [Step 4] Saving results...")
    
    # 4-1. 요약 데이터 생성
    summary_df = qc_executor.summarize(final_selection)
	print(summary_df)
    #if not summary_df.empty:
    #    # 4-2. Genomic Coordinate 복원 (Absolute Position)
    #    # Template상의 좌표 + Fetch 시작점(region_start_pos)
    #    genomic_starts = []
    #    genomic_ends = []
        
    #    for amp in final_selection:
    #        # Forward Start (5')
    #        g_start = region_start_pos + amp.forward.start
    #        # Reverse Start (5')
    #        g_end = region_start_pos + amp.reverse.start
            
    #        genomic_starts.append(g_start)
    #        genomic_ends.append(g_end)
            
    #    # 컬럼 삽입
    #    summary_df.insert(1, "Genomic_Fwd_Start", genomic_starts)
    #    summary_df.insert(2, "Genomic_Rev_Start", genomic_ends)
    #    summary_df.insert(0, "Chrom", args.chrom)

    # 4-3. 파일 쓰기
    if args.out.endswith(".xlsx"):
        summary_df.to_excel(args.out, index=False)
    else:
        summary_df.to_csv(args.out, index=False)

    print(f"✅ Done! Saved to: {args.out}")

if __name__ == "__main__":
    main()