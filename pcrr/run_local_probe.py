import sys
import os
import pandas as pd

# 1. src 폴더를 라이브러리 경로에 추가
sys.path.append(os.path.abspath("src"))

from src.pcr.config.loader import load_config
from src.pcr.designers.base import PrimerDesigner, ProbePrimerDesigner
from src.pcr.qc.executor import QCExecutor

def main():
    # -------------------------------------------------------------------------
    # [Step 1] Config 로드
    # -------------------------------------------------------------------------
    try:
        app_cfg = load_config() 
        print("✅ Config loaded successfully.")
    except Exception as e:
        print(f"❌ Config load failed: {e}")
        sys.exit(1)

    # -------------------------------------------------------------------------
    # [Step 2] 설정 (N수 증가 및 조건)
    # -------------------------------------------------------------------------
    overrides = {
        "PRIMER_PICK_INTERNAL_OLIGO": 1,
        
        # ✅ 후보군을 100개까지 늘려서 하나라도 건질 확률을 높임
        "PRIMER_NUM_RETURN": 100,
        
        # Probe 설정
        "max_probe_poly_g": 4,
        "max_probe_3_end_gc": 2,
        
        "PRIMER_INTERNAL_OPT_SIZE": 20, "PRIMER_INTERNAL_MIN_SIZE": 18, "PRIMER_INTERNAL_MAX_SIZE": 27,
        "PRIMER_INTERNAL_OPT_TM": 60.0, "PRIMER_INTERNAL_MIN_TM": 57.0, "PRIMER_INTERNAL_MAX_TM": 63.0,
        
        # Primer 설정
        "PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": 55.0, "PRIMER_MIN_TM": 50.0, "PRIMER_MAX_TM": 60.0,
        
        # 거리 제한 (너무 좁으면 디자인 실패 원인이 됨, 일단 3bp로 완화)
        "min_primer_probe_distance": 3,
    }

    primer_config = app_cfg.pcr_params.primer_kwargs
    probe_config = app_cfg.pcr_params.probe_kwargs

    # -------------------------------------------------------------------------
    # [Step 3] Primer Design Execution
    # -------------------------------------------------------------------------
    # ✅ [변경] 테스트를 위해 충분히 긴 서열 사용 (약 300bp)
    # 84bp로는 공간이 부족해 후보가 안 나옵니다.
    seq = (
        "CATGGACTCTCTGCACGAGTGCCCCCACATCCCCACCCCTGTGTGCATTGGGCACCGGGATGCACCCTCCTTCTCATCCCCACCACGCCAGGCTCCTGAGCCATCCCTCTTCTTCCAGGATCCCCCTGGAACTAGTATGGAGA"
    )
    print(seq)
    # Target 위치도 중간 지점으로 조정
    target_start = 110 
    target_len = 1
    use_probe = 1 

    print(f"🧪 Testing Sequence Length: {len(seq)} bp")
    print(f"🧪 Requesting {overrides['PRIMER_NUM_RETURN']} candidates...")

    if use_probe == 1:
        designer = ProbePrimerDesigner(
            seq, target_start, target_start + target_len, 
            config=primer_config, 
            probe_config=probe_config,
            overrides=overrides
        )
    else:
        designer = PrimerDesigner(
            seq, target_start, target_start + target_len, 
            config=primer_config, 
            overrides=overrides
        )

    try:
        results = designer.design()
        print(f"   -> Designed {len(results)} amplicons.")
    except Exception as e:
        print(f"❌ Design Failed: {e}")
        sys.exit(1)

    if not results:
        print("⚠️ No primers found. (Check constraints or sequence length)")
        sys.exit(0)

    # -------------------------------------------------------------------------
    # [Step 4] QC 실행
    # -------------------------------------------------------------------------
    print("\n🔍 [Step 3] Running QC Executor...")
    executor = QCExecutor(app_cfg.qc_params)
    checked_amplicons = executor.run_qc(results)

    # -------------------------------------------------------------------------
    # [Step 5] 검증 리포트
    # -------------------------------------------------------------------------
    print("\n📊 [Step 4] Validation Report")
    
    # PASS된 것들을 앞으로 정렬 (확인을 위해)
    checked_amplicons.sort(key=lambda x: x.is_qc_pass, reverse=True)
    
    df = executor.summarize(checked_amplicons)
    print(df.head(10).to_string(index=False)) # 너무 많으니 상위 10개만 요약 출력

    print(f"\n🔎 Detailed Verification (Top 3 Candidates):")
    
    # 상위 3개만 자세히 출력 (PASS된 것 우선)
    for i, amp in enumerate(checked_amplicons[:3]):
        status = "✅ PASS" if amp.is_qc_pass else "❌ FAIL"
        print(f"\n{'='*60}")
        print(f"[Candidate #{i+1}] {status} | Product Size: {amp.product_size}bp")
        print(f"{'='*60}")
        
        # 1. Structure & Map
        print(f"   [Structure & Map]")
        print(f"     Total Seq: {amp.sequence}")
        
        vis_map = ["."] * len(amp.sequence)
        
        # Fwd Map
        for k in range(len(amp.forward_primer.sequence)):
            if k < len(vis_map): vis_map[k] = ">"
            
        # Rev Map
        seq_len = len(amp.sequence)
        r_len = len(amp.reverse_primer.sequence)
        for k in range(seq_len - r_len, seq_len):
            if k >= 0: vis_map[k] = "<"
            
        # Probe Map
        if amp.probe:
            p_seq = amp.probe.sequence
            p_idx = amp.sequence.find(p_seq)
            if p_idx != -1:
                for k in range(p_idx, p_idx + len(p_seq)):
                     if k < len(vis_map): vis_map[k] = "P"

        print(f"     Ref Map  : {''.join(vis_map)}")
        print(f"                (>: Fwd, <: Rev, P: Probe, .: Template)")

        # 2. Sequences
        print(f"\n   [Sequences]")
        print(f"     Fwd: {amp.forward_primer.sequence}")
        if amp.probe: print(f"     Prb: {amp.probe.sequence}")
        print(f"     Rev: {amp.reverse_primer.sequence}")

        # 3. Tm Check
        p_tm = amp.probe.tm if amp.probe else 0.0
        f_tm = amp.forward_primer.tm
        r_tm = amp.reverse_primer.tm
        print(f"\n   [Tm Check]")
        print(f"     🔹 Probe Tm: {p_tm:.2f} C")
        print(f"     🔸 Fwd Tm:   {f_tm:.2f} C (Diff: {p_tm - f_tm:.2f})")
        print(f"     🔸 Rev Tm:   {r_tm:.2f} C (Diff: {p_tm - r_tm:.2f})")

        # 4. Thermo QC Details (★ Homodimer 포함 ★)
        t_data = amp.qc_status.thermo.data
        print(f"\n   [Thermo QC Details]")
        
        # Hairpin
        print(f"     🔹 Hairpin dG (Max: -4.0):")
        print(f"        - Fwd: {t_data.get('fwd_hairpin_dg', 0.0):.2f}")
        print(f"        - Rev: {t_data.get('rev_hairpin_dg', 0.0):.2f}")
        if 'probe_hairpin_dg' in t_data:
            print(f"        - Prb: {t_data.get('probe_hairpin_dg', 0.0):.2f}")

        # Homodimer (None Check 함수 사용)
        def fmt(val): return f"{val:.2f}" if val is not None else "0.00"
        
        f_homo = t_data.get('fwd_homodimer_dg')
        r_homo = t_data.get('rev_homodimer_dg')
        p_homo = t_data.get('probe_homodimer_dg')

        print(f"     🔹 Homodimer dG (Max: -9.0):")
        print(f"        - Fwd-Fwd: {fmt(f_homo)} {'❌ FAIL' if f_homo and f_homo < -9.0 else ''}")
        print(f"        - Rev-Rev: {fmt(r_homo)} {'❌ FAIL' if r_homo and r_homo < -9.0 else ''}")
        if p_homo is not None:
            print(f"        - Prb-Prb: {fmt(p_homo)} {'❌ FAIL' if p_homo < -9.0 else ''}")

        # Heterodimer
        print(f"     🔹 Heterodimer dG (Max: -9.0):")
        print(f"        - Fwd-Rev: {fmt(t_data.get('hetero_fr_dg'))} ({'OK' if t_data.get('hetero_fr_pass') else 'FAIL'})")
        
        if 'hetero_fp_dg' in t_data:
            print(f"        - Fwd-Prb: {fmt(t_data.get('hetero_fp_dg'))} ({'OK' if t_data.get('hetero_fp_pass') else 'FAIL'})")
        if 'hetero_rp_dg' in t_data:
            print(f"        - Rev-Prb: {fmt(t_data.get('hetero_rp_dg'))} ({'OK' if t_data.get('hetero_rp_pass') else 'FAIL'})")

        # 5. Specificity
        spec = amp.qc_status.specificity
        if spec.get('skipped'):
            print("\n   [Specificity] Skipped (Thermo Failed)")
        else:
            print(f"\n   [Specificity] {'OK' if spec['passed'] else 'FAIL'} (Off-targets: {spec.get('off_target_count')})")

if __name__ == "__main__":
    main()