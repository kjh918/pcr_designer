"""
scripts/run_pipeline.py
PCRFactory를 통한 통합 실행 스크립트.
"""
import argparse
import sys
import os
import pysam
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr import PCRFactory
from pcr.config.loader import ConfigLoader
from pcr.components.region import GenomicRegion


def parse_args():
    parser = argparse.ArgumentParser(description="PCR Pipeline")
    parser.add_argument("--assay",        type=str, required=True,
                        choices=["qpcr", "as_pcr", "ms_pcr"])
    parser.add_argument("--name",         type=str, required=True)
    parser.add_argument("--chrom",        type=str, required=True)
    parser.add_argument("--start",        type=int, required=True, help="1-based")
    parser.add_argument("--end",          type=int, required=True, help="1-based")
    parser.add_argument("--ref_genotype", type=str, default="")
    parser.add_argument("--alt_genotype", type=str, default="")
    parser.add_argument("--padding",      type=int, default=75)
    parser.add_argument("--ref",          type=str, default="hg38")
    parser.add_argument("--out",          type=str, default="results.xlsx")
    parser.add_argument("--top_k",        type=int, default=5)
    parser.add_argument("--num_search",   type=int, default=100)
    parser.add_argument("--no_qc",        action="store_true")
    return parser.parse_args()


def fetch_sequence(chrom, g_start, g_end, r_geno, a_geno, padding, fasta_path):
    with pysam.FastaFile(fasta_path) as fa:
        fetch_start = max(0, g_start - padding)
        fetch_end   = g_end + padding
        ref_seq     = fa.fetch(chrom, fetch_start, fetch_end).upper()
        rel_start   = g_start - fetch_start
        rel_end_ref = g_end   - fetch_start

        fetched = ref_seq[rel_start:rel_end_ref]
        if r_geno and fetched != r_geno:
            print(f"⚠️  Ref mismatch: genome='{fetched}' vs input='{r_geno}'")

        if a_geno:
            tmpl_seq  = ref_seq[:rel_start] + a_geno + ref_seq[rel_end_ref:]
            new_end   = rel_start + len(a_geno)
        else:
            tmpl_seq  = ref_seq
            new_end   = rel_end_ref

        return ref_seq, tmpl_seq, rel_start, new_end, fetch_start


def main():
    args = parse_args()

    config = ConfigLoader().load()
    if args.ref not in config.references:
        print(f"❌ Reference '{args.ref}' not found"); sys.exit(1)

    ref_path = config.references[args.ref].fasta_path

    print(f"🧬 [Step 1] Fetching sequence...")
    ref_seq, tmpl_seq, rel_start, rel_end, region_start = fetch_sequence(
        args.chrom, args.start - 1, args.end,
        args.ref_genotype, args.alt_genotype,
        args.padding, ref_path
    )

    # assay별 추가 kwargs
    extra = {}
    if args.assay == "as_pcr":
        extra = {"ref_genotype": args.ref_genotype, "alt_genotype": args.alt_genotype}

    print(f"🚀 [Step 2] Running {args.assay.upper()} pipeline...")
    factory = PCRFactory(config=config)
    output  = factory.run(
        assay_type=args.assay,
        name=args.name,
        template_sequence=tmpl_seq,
        reference_sequence=ref_seq,
        target_start=rel_start,
        target_end=rel_end,
        reference_name=args.ref,
        top_k=args.top_k,
        run_qc=not args.no_qc,
        overrides={"PRIMER_NUM_RETURN": args.num_search},
        **extra,
    )

    print(f"   Status : {output.status}")
    for msg in output.log_messages:
        print(f"   {msg}")

    if not output.amplicons:
        print("❌ No results."); sys.exit(1)

    # Genomic 좌표 주입
    for amp in output.amplicons:
        amp.region = GenomicRegion(
            chrom=args.chrom,
            start=region_start + amp.forward.start_index,
            end=region_start   + amp.reverse.start_index + len(amp.reverse.sequence),
            strand="+"
        )

    # 저장
    from pcr.qc.executor import QCExecutor
    qc  = QCExecutor(config)
    df  = qc.summarize(output.amplicons)

    if not df.empty:
        df.insert(0, "Chrom",       args.chrom)
        df.insert(1, "Ref_Allele",  args.ref_genotype)
        df.insert(2, "Alt_Allele",  args.alt_genotype)

    if args.out.endswith(".xlsx"):
        df.to_excel(args.out, index=False)
    else:
        df.to_csv(args.out, index=False)

    print(f"✅ Saved → {args.out}")


if __name__ == "__main__":
    main()