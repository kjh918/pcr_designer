# 명령어 맨 앞에 PYTHONPATH=. 을 추가해서 실행해 보세요!
#PYTHONPATH=. python3.11 scripts/design_qpcr.py \
#    --chrom chr7 \
#    --start 55174777 \
#    --end 55174777 \
#    --ref T \
#    --alt C \
#    --fasta /storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa \
#    --top_k 5
#PYTHONPATH=. python3.11 scripts/design_qpcr.py \
#    --chrom chr17 \
#    --start 39723640 \
#    --end 39723640 \
#    --ref G \
#    --alt A \
#    --genome hg38 \
#    --top_k 5 \
#    --base_config /tmp/tmppkkzed98.yaml

#PYTHONPATH=. python3.11 scripts/design_aspcr.py \
#    --chrom chr7 \
#    --start 55174777 \
#    --end 55174777 \
#    --ref T \
#    --alt C \
#    --fasta /storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa \
#    --top_k 5
#PYTHONPATH=. python3.11 scripts/design_mspcr.py \
#    --step 1 --seq "ATGCGCGCGATCGATCGATCGATCGCGC" --cpgs "5,26"
    
#PYTHONPATH=. python3.11 scripts/design_mspcr.py \
#    --step 2 --seq "GGTGCTTGGATCTGGCGCTTTTGGCACAGTCTACGAAGGT[CG]AGGGCCAGGTCCTGGGGTGGGCGGCCCCAGAGGATGGGGGCGGTGCCTGGAGGGGTGTGGTCGGCAGTTCTGATGGGAGGGGCAAGAGCTGGAGGCAGTGTTTGG" -k 5 --name "MGMT_Promoter"

PYTHONPATH=. python3.11 scripts/design_qc.py \
    --name test \
    --fwd GAAAATGACAAAGAACAGCTC \
    --rev TAGCACTTACCTGTA \
    #--probe CTGAAATCACTAAGCAGGAGA
