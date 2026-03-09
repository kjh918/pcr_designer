#!/usr/bin/env python3
import sys
import os
import json
import argparse
from typing import Dict, Any, List, Optional

# 상위 경로 추가 (pcr 패키지 인식용)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import load_pipeline_config
from pcr.factory import PCRFactory

# =====================================================================
# [Module 1] Bisulfite Conversion 시뮬레이션 (Step 1)
# =====================================================================
def simulate_bisulfite_conversion(raw_sequence: str, target_cpg_indices: Optional[List[int]] = None) -> Dict[str, Any]:
    """
    [Step 1 전용] 서열을 받아 M-Allele(Methylated)과 U-Allele(Unmethylated)로 변환합니다.
    - U-Allele: 모든 'C'가 'T'로 변환됨 (비메틸화 가정)
    - M-Allele: CpG 컨텍스트의 'C'는 'C'로 유지, 나머지 'C'는 'T'로 변환됨 (메틸화 가정)
    """
    seq = raw_sequence.strip().upper().replace(" ", "")
    
    u_allele = []
    m_allele = []
    cpg_positions = [] # 1-based index

    for i in range(len(seq)):
        base = seq[i]
        
        # U-Allele는 묻지도 따지지도 않고 C -> T 변환
        if base == 'C':
            u_allele.append('T')
        else:
            u_allele.append(base)
            
        # M-Allele는 CpG 여부 판단
        if base == 'C':
            if i < len(seq) - 1 and seq[i+1] == 'G':
                m_allele.append('C') # 메틸화되어 보호됨
                cpg_positions.append(i + 1)
            else:
                m_allele.append('T') # 메틸화 안 된 일반 C는 T로 변환됨
        else:
            m_allele.append(base)

    # 사용자가 특정 타겟 CpG를 지정하지 않았다면, 서열 내 전체 CpG를 타겟으로 간주
    valid_targets = target_cpg_indices if target_cpg_indices else cpg_positions

    return {
        "raw_sequence": seq,
        "m_allele": "".join(m_allele),
        "u_allele": "".join(u_allele),
        "total_cpg_count": len(cpg_positions),
        "target_cpg_count": len(valid_targets),
        "all_cpg_positions": cpg_positions,
        "target_cpg_positions": valid_targets
    }

# =====================================================================
# [Module 2] MS-PCR Primer Design 메인 파이프라인 (Step 2)
# =====================================================================
def design_manual_mspcr(
    design_name: str,
    raw_sequence: str,
    target_cpg_indices: Optional[List[int]] = None, 
    genome: str = "hg38",
    top_k: int = 5,
    base_yaml: str = "pcr/config/base_pcr.yaml",
    system_yaml: str = "pcr/config/system.yaml",
    qc_overrides: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    사용자가 입력한 Raw Sequence와 Target CpG를 기반으로 
    MS-PCR 프라이머(M-set, U-set)를 설계합니다.
    """
    assay_type = "mspcr"
    
    # 1. Module 1 호출하여 변환 데이터 확보
    conversion_info = simulate_bisulfite_conversion(raw_sequence, target_cpg_indices)
    
    # 타겟 구간 산출 (선택한 CpG들의 처음과 끝을 포함하는 범위)
    targets = conversion_info["target_cpg_positions"]
    print(targets)

    if not targets:
        return {"status": "fail", "reason": "No CpG sites found in the provided sequence."}
    
    # 0-based 구간 계산
    rel_start = min(targets) - 1
    rel_end = max(targets) + 1 # CG의 G까지 포함

    # 2. Config 로드
    user_overrides = {}
    if qc_overrides:
        user_overrides["qc_criteria"] = qc_overrides

    config = load_pipeline_config(base_yaml, system_yaml, assay_type, user_overrides=user_overrides)

    # 3. Factory 실행
    factory = PCRFactory(config)
    
    # MS-PCR 파이프라인으로 M과 U 두 개의 템플릿을 넘깁니다.
    templates = {
        "M": conversion_info["m_allele"],
        "U": conversion_info["u_allele"]
    }
    print(templates)
    output = factory.run(
        assay_type=assay_type,
        name=design_name,
        target_start=rel_start,
        target_end=rel_end,
        reference_name=genome,       
        template_sequence=conversion_info["raw_sequence"], # 원본
        templates=templates,         # 🔥 M/U 변환 서열 전달
        top_k=top_k,
        run_qc=True,                 
        overrides={"PRIMER_NUM_RETURN": 10},
        template_genomic_start=0,    
        target_strand="+"            
    )

    if output.status != "success":
        fail_reason = output.error_msg or "Failed to design MS-PCR primers."
        logs = getattr(output, "log_messages", [])
        return {
            "status": "fail", 
            "reason": fail_reason,
            "log_messages": "; ".join(logs) if isinstance(logs, list) else str(logs)
        }

    # 4. 결과 포맷팅 (Set 구조화)
    sets_dict = {}
    for amp in output.amplicons: 
        set_id = getattr(amp, "set_id", "UnknownSet")
        if set_id not in sets_dict:
            sets_dict[set_id] = []
        sets_dict[set_id].append(amp)

    results = []
    for rank, (set_id, amplicons_in_set) in enumerate(sets_dict.items(), start=1):
        set_qc_pass = all(getattr(amp, "is_qc_pass", False) for amp in amplicons_in_set)
        
        set_data = {
            "rank": rank,
            "set_id": set_id,
            "set_qc_pass": set_qc_pass,
            "alleles": {}
        }

        for amp in amplicons_in_set:
            allele_type = getattr(amp, "allele_type", "unknown") # "M" 또는 "U"
            alignment_data = getattr(amp, "alignment_visual", [])
            
            set_data["alleles"][allele_type] = {
                "id": amp.id,
                "metrics": {"pair_penalty": round(getattr(amp, "pair_penalty", 0), 3)},
                "qc_info": {"is_pass": amp.is_qc_pass, "fail_reason": getattr(amp, "qc_log", "None")},
                "amplicon_info": {
                    "sequence": amp.sequence if hasattr(amp, "sequence") else "-",
                    "tm": round(amp.tm, 2) if hasattr(amp, "tm") else 0, 
                    "gc": round(getattr(amp, "gc_percent", 0), 2)
                },
                "oligos": {
                    "forward": {
                        "sequence": amp.forward.sequence, 
                        "tm": round(amp.forward.tm, 2), 
                        "gc": round(getattr(amp.forward, "gc_percent", 0), 2),
                    },
                    "reverse": {
                        "sequence": amp.reverse.sequence, 
                        "tm": round(amp.reverse.tm, 2), 
                        "gc": round(getattr(amp.reverse, "gc_percent", 0), 2),
                    }
                },
                "alignment_text_block": "\n".join(alignment_data) if alignment_data else ""
            }
        
        results.append(set_data)

    return {
        "status": "success",
        "metadata": {
            "assay": "Manual-MS-PCR", 
            "design_name": design_name,
            "candidates_found": len(results),
            "genome_build_for_blast": genome
        },
        "conversion_info": conversion_info, 
        "target_info": {
            "template_length": len(conversion_info["raw_sequence"]),
            "target_range": f"{min(targets)}-{max(targets)}",
        },
        "results": results
    }

# =====================================================================
# CLI 실행부 (argparse 적용)
# =====================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MS-PCR Design Pipeline (Two-Step)")
    
    parser.add_argument("--step", type=int, choices=[1, 2], required=True, 
                        help="Step 1: Bisulfite Conversion only. Step 2: Full MS-PCR Design.")
    parser.add_argument("--name", type=str, default="MS-PCR_Manual_Target", help="Design Name / Project ID")
    parser.add_argument("--seq", type=str, required=True, help="Raw sequence (5' -> 3')")
    parser.add_argument("--cpgs", type=str, default="", help="Comma-separated 1-based indices of target CpGs (e.g., '24,58')")
    parser.add_argument("--genome", type=str, default="hg38", help="Reference genome for BLAST")
    parser.add_argument("-k", "--top_k", type=int, default=5, help="Number of candidate sets to return")
    parser.add_argument("--base_config", type=str, default="pcr/config/base_pcr.yaml")
    parser.add_argument("--system_config", type=str, default="pcr/config/system.yaml")

    args = parser.parse_args()

    # 타겟 CpG 인덱스 파싱
    target_cpg_list = []
    if args.cpgs:
        try:
            target_cpg_list = [int(x.strip()) for x in args.cpgs.split(",") if x.strip()]
        except ValueError:
            print(json.dumps({"status": "error", "reason": "Invalid --cpgs format. Please provide comma-separated integers."}))
            sys.exit(1)

    # 파이프라인 분기 처리
    if args.step == 1:
        # Step 1: Conversion 시뮬레이션 결과만 JSON으로 반환
        conversion_res = simulate_bisulfite_conversion(args.seq, target_cpg_list)
        print(json.dumps({"status": "success", "step": 1, "conversion_info": conversion_res}, indent=2))
        
    elif args.step == 2:
        # Step 2: 전체 프라이머 디자인 로직 수행
        design_res = design_manual_mspcr(
            design_name=args.name,
            raw_sequence=args.seq,
            target_cpg_indices=target_cpg_list,
            genome=args.genome,
            top_k=args.top_k,
            base_yaml=args.base_config,
            system_yaml=args.system_config
        )
        print(json.dumps(design_res, indent=2))