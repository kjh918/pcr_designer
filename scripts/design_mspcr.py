#!/usr/bin/env python3
import sys
import os
import json
import argparse
from typing import Dict, Any, List, Optional

# 상위 경로 추가
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import load_pipeline_config
from pcr.factory import PCRFactory
from pcr.designers.base.designer import BasePrimerDesigner 
from pcr.designers.ms_pcr.schema import MspcrDesignOutput

def simulate_bisulfite_conversion(clean_seq: str, target_indices: List[int]) -> Dict[str, Any]:
    """순수 서열과 타겟 위치를 받아 M/U Allele로 변환"""
    seq = clean_seq.upper()
    u_allele, m_allele, cpg_positions = [], [], []

    for i in range(len(seq)):
        base = seq[i]
        # U-Allele: 모든 C -> T
        u_allele.append('T' if base == 'C' else base)
        # M-Allele: CpG 'C'는 유지, 나머지는 'T'
        if base == 'C':
            if i < len(seq) - 1 and seq[i+1] == 'G':
                m_allele.append('C')
                cpg_positions.append(i + 1)
            else:
                m_allele.append('T')
        else:
            m_allele.append(base)

    return {
        "raw_sequence": seq,
        "m_allele": "".join(m_allele),
        "u_allele": "".join(u_allele),
        "total_cpg_count": len(cpg_positions),
        "target_cpg_count": len(target_indices),
        "all_cpg_positions": cpg_positions,
        "target_cpg_positions": target_indices
    }

def design_mspcr_primers(
    design_name: str,
    raw_sequence_with_brackets: str,
    genome: str = "hg38",
    top_k: int = 5,
    base_yaml: str = "pcr/config/base_pcr.yaml",
    system_yaml: str = "pcr/config/system.yaml",
    # 🔥 [수정] TypeError 해결을 위해 pcr_overrides 추가
    pcr_overrides: Optional[Dict[str, Any]] = None,
    qc_overrides: Optional[Dict[str, Any]] = None,
    window_size_3prime: int = 3,
    min_cpg_count: int = 1
) -> Dict[str, Any]:
    """
    MS-PCR 프라이머 설계 파이프라인.
    """
    assay_type = "mspcr"
    
    # 0. 대괄호 파싱
    clean_seq, ref_seq, target_indices = BasePrimerDesigner.parse_sequence_with_brackets(raw_sequence_with_brackets)
    
    if not target_indices:
        return {"status": "fail", "reason": "Target CpG [CG] not found in sequence."}

    # 1. Bisulfite 시뮬레이션
    conversion_info = simulate_bisulfite_conversion(clean_seq, target_indices)
    
    # 2. Config 로드 및 BLAST 스킵 설정
    user_overrides = {}
    if qc_overrides: user_overrides["qc_criteria"] = qc_overrides
    if pcr_overrides: user_overrides["pcr_params"] = pcr_overrides

    config = load_pipeline_config(base_yaml, system_yaml, assay_type, user_overrides=user_overrides)
    
    # BLAST 스킵 로직 (reference_name이 'none'인 경우)
    if hasattr(config, "system") and hasattr(config.system, "paths"):
        if str(genome).lower() == "none":
            config.system.paths.blast_db_path = None

    # 3. PCRFactory 실행
    factory = PCRFactory(config)
    
    templates = {
        "M": conversion_info["m_allele"],
        "U": conversion_info["u_allele"]
    }
    
    # Extra Input Kwargs (Designer 내부로 전달됨)
    extra_input_kwargs = {
        "templates": templates,
        "target_cpg_indices": target_indices,
        "window_size_3prime": window_size_3prime,
        "min_cpg_count": min_cpg_count
    }

    # Factory Run
    output = factory.run(
        assay_type=assay_type,
        name=design_name,
        target_start=min(target_indices) - 1,
        target_end=max(target_indices),
        reference_name=genome,       
        template_sequence=conversion_info["raw_sequence"], 
        top_k=top_k,
        run_qc=True,
        **extra_input_kwargs
    )

    if output.status != "success":
        return {"status": "fail", "reason": output.error_msg or "MS-PCR design failed."}

    # 4. 결과 포맷팅 (MspcrDesignOutput 스키마 사용)
    mspcr_output = MspcrDesignOutput(
        amplicons=output.amplicons,
        status="success",
        metadata={
            "assay": "MS-PCR", 
            "genome": genome, 
            "project": design_name,
            "reference_genome": genome
        }
    )

    result_dict = mspcr_output.to_frontend_dict()
    
    # UI 추가 정보 병합
    result_dict["conversion_info"] = conversion_info
    result_dict["target_info"] = {
        "template_length": len(clean_seq),
        "target_range": f"{min(target_indices)}-{max(target_indices)}",
    }

    return result_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", default="MS-PCR_CLI")
    parser.add_argument("--seq", required=True)
    parser.add_argument("--genome", default="hg38")
    args = parser.parse_args()
    print(json.dumps(design_mspcr_primers(args.name, args.seq, args.genome), indent=2))