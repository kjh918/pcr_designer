# 명령어 맨 앞에 PYTHONPATH=. 을 추가해서 실행해 보세요!
PYTHONPATH=. python scripts/design_qpcr.py \
    --chrom chr17 \
    --start 39723627 \
    --end 39723628 \
    --fasta /Users/kimjihoon/Documents/reference_genome/hg38.fa \
    --top_k 5