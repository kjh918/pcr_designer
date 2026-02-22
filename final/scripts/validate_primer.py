#!/usr/bin/env python3
"""
scripts/validate_primers.py
기존 프라이머/프로브 세트를 검증(In-Silico PCR 및 열역학 QC)하는 독립 실행형 스크립트.
다른 모듈에서 `validate_primers` 함수를 직접 import 하여 사용할 수 있습니다.
"""
import sys
import os
import json
import argparse
from typing import Dict, Any, Optional

# 프로젝트 최상단 디렉토리를 경로에 추가하여 pcr 모듈 인식
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import load_pipeline_config
from pcr.components.amplicon import Amplicon
from pcr.components.primer import Primer, Probe
from pcr.qc.blast import BlastSpecificityChecker
from pcr.qc.thermo import ThermoChecker

def validate_primers(
    fwd_seq: str, 
    rev_seq: str, 
    probe_seq: Optional[str] = None, 
    assay_type: str = "qpcr",
    yaml_path: str = "pcr/config/base_pcr.yaml"
) -> Dict[str, Any]:
    """
    [MODIFIED] Class 래퍼를 제거하고 독립적인 함수로 분리했습니다.
    입력된 서열을 바탕으로 BLAST(특이성)와 Thermo(열역학) 검증을 수행합니다.
    """
    fwd_seq = fwd_seq.upper()
    rev_seq = rev_seq.upper()
    probe_seq = probe_seq.upper() if probe_seq else None

    # 1. 설정 로드 및 도구 초기화
    try:
        config = load_pipeline_config(yaml_path, assay_type=assay_type)
        blast_checker = BlastSpecificityChecker(config)
        thermo_checker = ThermoChecker(config.qc_criteria)
    except Exception as e:
        return {"status": "error", "reason": f"Failed to load config or initialize checkers: {e}"}

    # ---------------------------------------------------------
    # Step 1: BLAST (In-Silico PCR)
    # ---------------------------------------------------------
    hits = blast_checker._run_blast_set(fwd_seq, rev_seq, probe_seq)
    
    f_hits = [h for h in hits if h.qseqid == "Fwd"]
    r_hits = [h for h in hits if h.qseqid == "Rev"]
    p_hits = [h for h in hits if h.qseqid == "Probe"]

    found_amplicons = blast_checker._find_amplicons(f_hits, r_hits)
    
    if not found_amplicons:
        return {"status": "fail", "reason": "No valid amplification product found via BLAST.", "amplicons": []}

    # ---------------------------------------------------------
    # Step 2: Amplicon 객체화 및 Thermo QC
    # ---------------------------------------------------------
    results = []
    for idx, cand in enumerate(found_amplicons):
        # 도메인 객체 생성
        amp = Amplicon(
            id=f"Valid_{idx}",
            forward=Primer(sequence=fwd_seq),
            reverse=Primer(sequence=rev_seq),
            probe=Probe(sequence=probe_seq) if probe_seq else None,
            reference_id=cand.chrom,
            target_start_index=cand.start,
            target_end_index=cand.end,
            template_sequence="" # 필요시 fetch
        )
        
        # Probe 결합 위치 확인
        probe_binds = False
        if probe_seq:
            probe_binds = blast_checker._check_probe_binding(cand, p_hits)
            
        # Thermo QC 계산
        thermo_result = thermo_checker._check_single(amp)
        
        results.append({
            "chrom": cand.chrom,
            "start": cand.start,
            "end": cand.end,
            "product_size": cand.product_size,
            "probe_binds": probe_binds,
            "thermo_passed": thermo_result["passed"],
            "thermo_details": thermo_result["data"],
            "fail_reason": thermo_result["fail_reason"]
        })

    # ---------------------------------------------------------
    # Step 3: 최종 결과 포맷팅
    # ---------------------------------------------------------
    return {
        "status": "success",
        "assay_type": assay_type,
        "blast_summary": {
            "total_amplicons_found": len(found_amplicons),
            "intended_target_found": len(found_amplicons) == 1,
            "off_target_count": max(0, len(found_amplicons) - 1)
        },
        "amplicons": results
    }

# =====================================================================
# CLI 실행 모드 (터미널에서 직접 실행할 때 작동)
# =====================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Existing Primer/Probe Validation Tool")
    parser.add_argument("-f", "--fwd", required=True, help="Forward Primer Sequence")
    parser.add_argument("-r", "--rev", required=True, help="Reverse Primer Sequence")
    parser.add_argument("-p", "--probe", default=None, help="Probe Sequence (Optional)")
    parser.add_argument("-a", "--assay", default="qpcr", help="Assay Type (e.g., qpcr, ms_pcr, base)")
    parser.add_argument("-c", "--config", default="../pcr/config/base_pcr.yaml", help="Path to base_pcr.yaml")
    
    args = parser.parse_args()

    # 함수 직접 호출
    result = validate_primers(
        fwd_seq=args.fwd,
        rev_seq=args.rev,
        probe_seq=args.probe,
        assay_type=args.assay,
        yaml_path=args.config
    )

    # 결과를 JSON 형태로 예쁘게 출력 (다른 시스템에서 파이프(|)로 넘겨받기 좋음)
    print(json.dumps(result, indent=2, ensure_ascii=False))