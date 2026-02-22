#!/usr/bin/env python3
"""
scripts/design_qpcr.py
좌표(chrom, start, end)를 입력받아 TaqMan qPCR 프라이머/프로브를 설계하고,
QC(열역학/특이성) 및 랭킹이 적용된 최종 후보군을 반환하는 독립 실행형 스크립트입니다.
(최신 PCRFactory 아키텍처 적용)
"""
import sys
import os
import json
import argparse
from typing import Dict, Any

try:
    import pysam
except ImportError:
    pysam = None

# 프로젝트 최상단 디렉토리를 경로에 추가
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import load_pipeline_config
from pcr.factory import PCRFactory

def fetch_template_sequence(chrom: str, start: int, end: int, padding: int, fasta_path: str):
    """Reference FASTA에서 패딩이 포함된 Template 서열을 추출하고 상대 좌표를 계산합니다."""
    if not pysam:
        raise ImportError("pysam 라이브러리가 설치되어 있지 않습니다. (pip install pysam)")
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Reference FASTA 파일을 찾을 수 없습니다: {fasta_path}")

    fetch_start = max(0, start - padding)
    fetch_end = end + padding

    with pysam.FastaFile(fasta_path) as fasta:
        try:
            seq = fasta.fetch(chrom, fetch_start, fetch_end).upper()
        except KeyError:
            # chr 접두사 처리
            alt_chrom = chrom.replace("chr", "") if "chr" in chrom else f"chr{chrom}"
            seq = fasta.fetch(alt_chrom, fetch_start, fetch_end).upper()

    # 추출된 서열 내에서의 실제 타겟 상대 좌표
    relative_target_start = start - fetch_start
    relative_target_end = relative_target_start + (end - start)

    return seq, relative_target_start, relative_target_end

def design_qpcr_primers(
    chrom: str,
    start: int,
    end: int,
    fasta_path: str,
    padding: int = 150,
    top_k: int = 5,
    base_yaml: str = "pcr/config/base_pcr.yaml",
    system_yaml: str = "pcr/config/system.yaml"
) -> Dict[str, Any]:
    
    assay_type = "qpcr"
    task_name = f"{chrom}_{start}_{end}"

    # 1. Config 로드 (base와 system 분리 병합)
    try:
        config = load_pipeline_config(
            base_yaml_path=base_yaml,
            system_yaml_path=system_yaml,
            assay_type=assay_type
        )
    except Exception as e:
        return {"status": "error", "reason": f"Config load failed: {e}"}

    # 2. Reference 서열 추출
    try:
        template_seq, rel_start, rel_end = fetch_template_sequence(
            chrom, start, end, padding, fasta_path
        )
        print(template_seq, rel_start, rel_end)
    except Exception as e:
        return {"status": "error", "reason": f"Sequence extraction failed: {e}"}

    # ---------------------------------------------------------
    # 3. Factory 단일 실행 (Design -> QC -> Rank 원패스)
    # ---------------------------------------------------------
    factory = PCRFactory(config)
    
    try:
        output = factory.run(
            assay_type=assay_type,
            name=task_name,
            target_start=rel_start,
            target_end=rel_end,
            reference_name=os.path.basename(fasta_path),
            template_sequence=template_seq,
            top_k=top_k,
            run_qc=True,  # ★ 여기서 True를 주면 내부에서 BLAST, Thermo QC가 전부 실행됨
            overrides={"PRIMER_NUM_RETURN": 100} # Primer3 초기 생성 수
        )
    except Exception as e:
         return {"status": "error", "reason": f"Factory execution failed: {e}"}

    # 실패 처리
    if output.status != "success" or not output.amplicons:
        return {
            "status": "fail", 
            "reason": output.error_msg or "Design or QC failed.", 
            "log": output.log_messages
        }

    # ---------------------------------------------------------
    # 4. 결과 포맷팅
    # ---------------------------------------------------------
    results = []
    for amp in output.amplicons:
        results.append({
            "id": amp.id,
            "total_penalty": round(getattr(amp, 'total_penalty', amp.pair_penalty), 3),
            "forward": {"seq": amp.forward.sequence, "tm": round(amp.forward.tm, 2)},
            "reverse": {"seq": amp.reverse.sequence, "tm": round(amp.reverse.tm, 2)},
            "probe": {"seq": amp.probe.sequence, "tm": round(amp.probe.tm, 2)} if amp.probe else None,
            "qc_passed": amp.is_qc_pass,
            "thermo_stats": getattr(amp, 'thermo_stats', {}),
            "blast_stats": getattr(amp, 'blast_stats', {})
        })

    return {
        "status": "success",
        "task_name": task_name,
        "assay_type": assay_type,
        "log": output.log_messages,
        "summary": {
            "final_selected": len(output.amplicons)
        },
        "amplicons": results
    }

# =====================================================================
# CLI 실행 모드
# =====================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="qPCR Primer Design Tool (End-to-End)")
    parser.add_argument("-c", "--chrom", required=True, help="Chromosome (e.g., chr1)")
    parser.add_argument("-s", "--start", type=int, required=True, help="Target start coordinate")
    parser.add_argument("-e", "--end", type=int, required=True, help="Target end coordinate")
    parser.add_argument("-f", "--fasta", required=True, help="Path to Reference FASTA file")
    parser.add_argument("-p", "--padding", type=int, default=75, help="Padding around target (bp)")
    parser.add_argument("-k", "--top_k", type=int, default=5, help="Number of final candidates to return")
    parser.add_argument("--base_config", default="/Users/kimjihoon/Downloads/GdriveBackup/Projects/pcr_designer/final/pcr/config/base_pcr.yaml", help="Path to base_pcr.yaml")
    parser.add_argument("--system_config", default="/Users/kimjihoon/Downloads/GdriveBackup/Projects/pcr_designer/final/pcr/config/system.yaml", help="Path to system.yaml")
    
    args = parser.parse_args()

    result = design_qpcr_primers(
        chrom=args.chrom,
        start=args.start,
        end=args.end,
        fasta_path=args.fasta,
        padding=args.padding,
        top_k=args.top_k,
        base_yaml=args.base_config,
        system_yaml=args.system_config
    )

    print(json.dumps(result, indent=2, ensure_ascii=False))