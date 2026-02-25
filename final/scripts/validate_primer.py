#!/usr/bin/env python3
import sys
import os
import json
from typing import Dict, Any, List, Optional

# 상위 경로 추가 (pcr 패키지 인식용)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import load_pipeline_config
from pcr.factory import PCRFactory

def design_manual_qpcr(
    design_name: str,
    raw_sequence: str,
    target_start: int,  # 1-based start
    target_end: int,    # 1-based end
    genome: str = "hg38",
    top_k: int = 5,
    base_yaml: str = "pcr/config/base_pcr.yaml",
    system_yaml: str = "pcr/config/system.yaml",
    qc_overrides: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    사용자가 입력한 Raw Sequence 내에서 특정 Target 구간을 포함하는 
    qPCR 프라이머 및 프로브를 설계합니다.
    """
    
    assay_type = "qpcr"
    # 서열 정제 (공백 및 줄바꿈 제거)
    template = raw_sequence.strip().upper().replace(" ", "")
    
    # 1. Config 로드
    user_overrides = {}
    if qc_overrides:
        user_overrides["qc_criteria"] = qc_overrides

    config = load_pipeline_config(base_yaml, system_yaml, assay_type, user_overrides=user_overrides)

    # 2. 0-based 좌표 변환
    # 사용자가 150-155를 넣었다면 인덱스로는 149:155가 됨
    rel_start = target_start - 1
    rel_end = target_end
    print(1)
    # 3. Factory 실행
    factory = PCRFactory(config)
    
    # 매뉴얼 설계에서는 reference_sequence와 template_sequence가 동일하거나
    # 변이 전/후 서열을 직접 생성해야 함 (여기서는 입력 서열을 타겟으로 함)
    output = factory.run(
        assay_type=assay_type,
        name=design_name,
        target_start=rel_start,
        target_end=rel_end,
        reference_name=genome,       # BLAST 수행을 위해 필요
        template_sequence=template,  # 설계 대상 서열
        reference_sequence=template, # 비교 대상 서열
        top_k=top_k,
        run_qc=True,                 # QC 및 BLAST 활성화
        overrides={"PRIMER_NUM_RETURN": 10},
        template_genomic_start=0,    # 매뉴얼 서열이므로 0부터 시작으로 간주
        target_strand="+"            # 매뉴얼 입력은 항상 Plus strand 기준
    )
    print(output)
    if output.status != "success":
        fail_reason = output.error_msg or "Failed to design primers/probes."
        logs = getattr(output, "log_messages", [])
        return {
            "status": "fail", 
            "reason": fail_reason,
            "log_messages": "; ".join(logs) if isinstance(logs, list) else str(logs)
        }

    # 4. 결과 포맷팅 (기존 코드와 동일 구조)
    def get_position_info(seq: str, is_reverse: bool = False):
        if not seq: return None
        # DNA 서열 내에서 위치 탐색
        search_seq = seq.upper()
        # 역방향 프라이머의 경우 템플릿 서열 내에서 상보 서열 위치를 찾아야 함
        if is_reverse:
            # 역상보 변환 (간이 함수)
            rev_comp = search_seq.translate(str.maketrans('ATGC', 'TACG'))[::-1]
            rel_idx = template.find(rev_comp)
        else:
            rel_idx = template.find(search_seq)

        if rel_idx == -1:
            return {"index_start": None, "index_end": None, "strand": "?"}
        
        return {
            "index_start": rel_idx + 1,
            "index_end": rel_idx + len(seq),
            "strand": "+" if not is_reverse else "-"
        }

    results = []
    for rank, amp in enumerate(output.amplicons, start=1):
        alignment_data = getattr(amp, "alignment_visual", ["Alignment data not available."])
        
        results.append({
            "rank": rank,
            "id": amp.id,
            "metrics": {"pair_penalty": round(amp.pair_penalty, 3)},
            "qc_info": {"is_pass": amp.is_qc_pass, "fail_reason": getattr(amp, "qc_fail_reason", "None")},
            "amplicon_info": {
                "sequence": amp.sequence, 
                "length": amp.product_size,
                "tm": round(amp.tm, 2), 
                "gc": amp.gc
            },
            "oligos": {
                "forward": {
                    "sequence": amp.forward.sequence, 
                    "length": len(amp.forward.sequence),
                    "tm": round(amp.forward.tm, 2), 
                    "gc": round((amp.forward.sequence.count('G') + amp.forward.sequence.count('C')) / len(amp.forward.sequence) * 100, 2),
                    **get_position_info(amp.forward.sequence, False)
                },
                "reverse": {
                    "sequence": amp.reverse.sequence, 
                    "length": len(amp.reverse.sequence),
                    "tm": round(amp.reverse.tm, 2), 
                    "gc": round((amp.reverse.sequence.count('G') + amp.reverse.sequence.count('C')) / len(amp.reverse.sequence) * 100, 2),
                    **get_position_info(amp.reverse.sequence, True)
                },
                "probe": {
                    "sequence": amp.probe.sequence, 
                    "length": len(amp.probe.sequence),
                    "tm": round(amp.probe.tm, 2), 
                    "gc": round((amp.probe.sequence.count('G') + amp.probe.sequence.count('C')) / len(amp.probe.sequence) * 100, 2),
                    **get_position_info(amp.probe.sequence, False) # Probe는 일반적으로 Forward 방향
                } if amp.probe else None
            },
            "alignment_text_block": "\n".join(alignment_data)
        })

    return {
        "status": "success",
        "metadata": {
            "assay": "Manual-TaqMan-qPCR", 
            "design_name": design_name,
            "candidates_found": len(results),
            "genome_build_for_blast": genome
        },
        "target_info": {
            "template_length": len(template),
            "target_range": f"{target_start}-{target_end}",
            "template_sequence": template 
        },
        "results": results
    }

# 실행 예시
if __name__ == "__main__":
    test_seq = "ATGC..." # 실제 긴 서열 입력
    res = design_manual_qpcr(
        design_name="Test_Manual_Project",
        raw_sequence=test_seq,
        target_start=150,
        target_end=155
    )
    print(json.dumps(res, indent=2))