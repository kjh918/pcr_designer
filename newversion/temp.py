import json
from pcr.components import GenomicRegion, Primer, Probe, Amplicon

def run_test():
    print("=== PCR Components Test Start ===\n")

    # ------------------------------------------------------------------
    # 1. 가상 데이터 준비 (Scenario Setup)
    # ------------------------------------------------------------------
    # 길이 100bp 서열
    # Index 20: C (Ref) -> T (Template) : 의도치 않은 변이 (Mismatch)
    # Index 50: G (Ref) -> A (Template) : 타겟 변이 (Target SNP)
    
    # 기본 골격 (A로 채움)
    base_seq = ["A"] * 100
    
    # Reference 생성
    ref_list = base_seq.copy()
    ref_list[20] = "C"
    ref_list[50] = "G" # Wildtype
    reference_seq = "".join(ref_list)

    # Template (Target) 생성
    tmpl_list = base_seq.copy()
    tmpl_list[20] = "T" # Mismatch 발생 지점
    tmpl_list[50] = "A" # Target SNP 지점
    template_seq = "".join(tmpl_list)

    # Genomic Region (가상 좌표: chr1:1000-1100)
    template_region = GenomicRegion(chrom="chr1", start=1000, end=1100, strand="+", sequence=template_seq)

    print(f"Reference: {reference_seq[:60]}...")
    print(f"Template : {template_seq[:60]}...")
    print("-" * 60)

    # ------------------------------------------------------------------
    # 2. Primer3 결과 Mocking (가짜 결과 데이터)
    # ------------------------------------------------------------------
    # Forward: 0~20bp / Reverse: 80~100bp / Probe: 45~55bp (SNP 포함)
    primer3_result = {
        # Forward Primer (Left)
        "PRIMER_LEFT_0_SEQUENCE": "AAAAAAAAAAAAAAAAAAAA", # index 0, len 20
        "PRIMER_LEFT_0": [0, 20], 
        "PRIMER_LEFT_0_TM": 55.5,
        "PRIMER_LEFT_0_GC_PERCENT": 0.0,
        "PRIMER_LEFT_0_PENALTY": 0.1,
        "PRIMER_LEFT_0_HAIRPIN_TH": 0.0,
        "PRIMER_LEFT_0_SELF_ANY_TH": 0.0,

        # Reverse Primer (Right) - Primer3는 Right의 경우 3' end index를 줌
        "PRIMER_RIGHT_0_SEQUENCE": "TTTTTTTTTTTTTTTTTTTT", # index 80, len 20
        "PRIMER_RIGHT_0": [99, 20], # 0-based index 99에서 뒤로 20bp
        "PRIMER_RIGHT_0_TM": 55.0,
        "PRIMER_RIGHT_0_GC_PERCENT": 0.0,
        "PRIMER_RIGHT_0_PENALTY": 0.2,
        "PRIMER_RIGHT_0_HAIRPIN_TH": 0.0,
        "PRIMER_RIGHT_0_SELF_ANY_TH": 0.0,

        # Probe (Internal)
        "PRIMER_INTERNAL_0_SEQUENCE": "AAAAAGAAAA", # index 45, len 10
        "PRIMER_INTERNAL_0": [45, 10],
        "PRIMER_INTERNAL_0_TM": 60.0,
        "PRIMER_INTERNAL_0_GC_PERCENT": 10.0,
        "PRIMER_INTERNAL_0_PENALTY": 0.05,
        "PRIMER_INTERNAL_0_HAIRPIN_TH": 0.0,
        "PRIMER_INTERNAL_0_SELF_ANY_TH": 0.0,
    }

    # ------------------------------------------------------------------
    # 3. Component 생성 (Factory Method 사용)
    # ------------------------------------------------------------------
    # Forward Primer 생성
    fwd = Primer.from_primer3(primer3_result, rank=0, role_key="LEFT", template_region=template_region)
    
    # Reverse Primer 생성
    rev = Primer.from_primer3(primer3_result, rank=0, role_key="RIGHT", template_region=template_region)
    
    # Probe 생성 (Target Metadata 주입)
    # Target SNP 위치는 Template 기준 index 50이라고 가정
    probe = Probe.from_primer3(
        primer3_result, rank=0, template_region=template_region,
        target_id="rs12345",
        target_type="SNP",
        target_status="Alt",
        population="EAS" # 추가 메타데이터 테스트
    )

    if not fwd or not rev or not probe:
        print("Error: Failed to create primer objects.")
        return

    # ------------------------------------------------------------------
    # 4. Amplicon 조립 및 분석 수행
    # ------------------------------------------------------------------
    print("Building Amplicon and Analyzing Variations...")
    
    amplicon = Amplicon(
        forward=fwd,
        reverse=rev,
        probe=probe,
        template_sequence=template_seq,
        # 핵심: 타겟 영역 지정 (index 50~51)
        target_start_index=50,
        target_end_index=51,
        reference_sequence=reference_seq,
        reference_id="hg38",
        is_qc_pass=True
    )

    # ------------------------------------------------------------------
    # 5. 결과 확인 (JSON Output)
    # ------------------------------------------------------------------
    result_dict = amplicon.to_dict()
    
    # 보기 좋게 출력
    print("\n=== Final Result JSON ===")
    print(json.dumps(result_dict, indent=4))

    # ------------------------------------------------------------------
    # 6. 검증 포인트 확인
    # ------------------------------------------------------------------
    print("\n=== Validation Checks ===")
    
    # Check 1: Mismatch Count (Index 20의 'C>T' 하나만 잡혀야 함)
    mismatch = result_dict['mismatch_count']
    print(f"1. Unintended Mismatch Count: {mismatch} (Expected: 1)")
    
    # Check 2: Target Change (Index 50의 'G>A'는 타겟으로 잡혀야 함)
    target_cnt = result_dict['target_change_count']
    print(f"2. Target Change Count      : {target_cnt} (Expected: 1)")
    
    # Check 3: Changes Detail
    changes = result_dict['changes']
    print("3. Changes Detail:")
    for c in changes:
        status = "✅ TARGET" if c['on_target'] else "❌ MISMATCH"
        print(f"   - [{status}] Pos: {c['pos']}, {c['change']}, Conversion: {c['is_conversion']}")

if __name__ == "__main__":
    run_test()