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
# 🔥 대괄호 파싱을 위해 BasePrimerDesigner 호출
from pcr.designers.base.designer import BasePrimerDesigner 

# =====================================================================
# [Module 1] Bisulfite Conversion 시뮬레이션 (Step 1)
# =====================================================================
def simulate_bisulfite_conversion(clean_seq: str, target_indices: List[int]) -> Dict[str, Any]:
    """
    [Step 1 전용] 순수 서열과 추출된 타겟 위치를 받아 M/U Allele로 변환합니다.
    - U-Allele: 모든 'C'가 'T'로 변환됨 (비메틸화 가정)
    - M-Allele: CpG 컨텍스트의 'C'는 'C'로 유지, 나머지 'C'는 'T'로 변환됨 (메틸화 가정)
    """
    seq = clean_seq.upper()
    
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

    return {
        "raw_sequence": seq,
        "m_allele": "".join(m_allele),
        "u_allele": "".join(u_allele),
        "total_cpg_count": len(cpg_positions),
        "target_cpg_count": len(target_indices),
        "all_cpg_positions": cpg_positions,
        "target_cpg_positions": target_indices # 대괄호로 지정된 진짜 타겟 위치
    }

# =====================================================================
# [Module 2] MS-PCR Primer Design 메인 파이프라인 (Step 2)
# =====================================================================
def design_mspcr_primers(
    design_name: str,
    raw_sequence_with_brackets: str, # 🔥 숫자 인덱스 대신 대괄호 포함 서열
    genome: str = "hg38",
    top_k: int = 5,
    window_size_3prime: int = 4,     # 🔥 [NEW] 3' 말단 윈도우 사이즈 파라미터 추가
    min_cpg_count: int = 1,          # 🔥 [NEW] 프라이머 내 최소 CpG 개수 파라미터 추가
    base_yaml: str = "pcr/config/base_pcr.yaml",
    system_yaml: str = "pcr/config/system.yaml",
    qc_overrides: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    사용자가 입력한 서열(대괄호 타겟 지정)을 기반으로 
    MS-PCR 프라이머(M-set, U-set)를 설계합니다.
    """
    assay_type = "mspcr"
    
    # 0. 대괄호 파싱을 통해 순수 서열과 타겟 위치 자동 추출
    clean_seq, target_indices = BasePrimerDesigner.parse_sequence_with_brackets(raw_sequence_with_brackets)
    
    if not target_indices:
        return {"status": "fail", "reason": "Target not found. Please wrap your target CpG with brackets, e.g., [CG]"}

    # 1. Module 1 호출하여 변환 데이터(M/U 유무) 확보
    conversion_info = simulate_bisulfite_conversion(clean_seq, target_indices)
    
    # 0-based 구간 계산
    rel_start = min(target_indices) - 1
    rel_end = max(target_indices)

    # 2. Config 로드
    user_overrides = {}
    if qc_overrides:
        user_overrides["qc_criteria"] = qc_overrides

    config = load_pipeline_config(base_yaml, system_yaml, assay_type, user_overrides=user_overrides)

    # 3. Factory 실행
    factory = PCRFactory(config)
    
    templates = {
        "M": conversion_info["m_allele"],
        "U": conversion_info["u_allele"]
    }
    
    # Pydantic Input으로 넘겨주기 위한 추가 인자 (MS-PCR 전용)
    extra_input_kwargs = {
        "target_cpg_indices": target_indices,
        "templates": templates,
        "window_size_3prime": window_size_3prime, # 🔥 schema.py로 전달
        "min_cpg_count": min_cpg_count            # 🔥 schema.py로 전달
    }

    output = factory.run(
        assay_type=assay_type,
        name=design_name,
        target_start=rel_start,
        target_end=rel_end,
        reference_name=genome,       
        template_sequence=conversion_info["raw_sequence"], # 원본
        top_k=top_k,
        run_qc=True,                 
        overrides={"PRIMER_NUM_RETURN": 20},
        template_genomic_start=0,    
        target_strand="+",
        **extra_input_kwargs         # 🔥 M/U 변환 서열 및 신규 파라미터 주입
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
            "target_range": f"{min(target_indices)}-{max(target_indices)}",
        },
        "results": results
    }

# =====================================================================
# CLI 실행부 (argparse 적용)
# =====================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MS-PCR Design Pipeline (Two-Step with Brackets)")
    
    parser.add_argument("--step", type=int, choices=[1, 2], required=True, 
                        help="Step 1: Bisulfite Conversion only. Step 2: Full MS-PCR Design.")
    parser.add_argument("--name", type=str, default="MS-PCR_Manual_Target", help="Design Name / Project ID")
    parser.add_argument("--seq", type=str, required=True, help="Sequence with brackets, e.g., ATGC[CG]ATGC")
    parser.add_argument("--genome", type=str, default="hg38", help="Reference genome for BLAST")
    parser.add_argument("-k", "--top_k", type=int, default=10, help="Number of candidate sets to return")
    parser.add_argument("--window", type=int, default=4, help="Allowed window size at the 3' end") # 🔥 CLI 인자 추가
    parser.add_argument("--min_cpg", type=int, default=1, help="Minimum number of CpG sites required") # 🔥 CLI 인자 추가
    parser.add_argument("--base_config", type=str, default="pcr/config/base_pcr.yaml")
    parser.add_argument("--system_config", type=str, default="pcr/config/system.yaml")

    args = parser.parse_args()

    if args.step == 1:
        # Step 1: 파싱 후 Conversion 시뮬레이션 결과만 반환
        clean_seq, target_indices = BasePrimerDesigner.parse_sequence_with_brackets(args.seq)
        if not target_indices:
            print(json.dumps({"status": "fail", "reason": "Missing brackets [] in sequence."}, indent=2))
            sys.exit(1)
            
        conversion_res = simulate_bisulfite_conversion(clean_seq, target_indices)
        print(json.dumps({"status": "success", "step": 1, "conversion_info": conversion_res}, indent=2))
        
    elif args.step == 2:
        # Step 2: 전체 프라이머 디자인
        design_res = design_mspcr_primers(
            design_name=args.name,
            raw_sequence_with_brackets=args.seq,
            genome=args.genome,
            top_k=args.top_k,
            window_size_3prime=args.window, # 🔥 전달
            min_cpg_count=args.min_cpg,     # 🔥 전달
            base_yaml=args.base_config,
            system_yaml=args.system_config
        )
        print(json.dumps(design_res, indent=2))