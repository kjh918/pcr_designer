import sys
import os
import pandas as pd

# 1. src 폴더를 라이브러리 경로에 추가
sys.path.append(os.path.abspath("src"))

from src.pcr.config.loader import load_config
from src.pcr.designers.base import BasePrimerDesigner, ProbePrimerDesigner
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
	# [Step 2] Probe Design을 위한 설정 (Overrides)
	# -------------------------------------------------------------------------
	# Primer3에게 "Probe도 뽑아줘"라고 명령하려면 PRIMER_PICK_INTERNAL_OLIGO = 1 이어야 합니다.
	overrides = {
		# 1. Probe 활성화 (필수)
		"PRIMER_PICK_INTERNAL_OLIGO": 1,
		
		# 2. Probe 조건 설정 (Primer와 다르게 설정 가능)
		"PRIMER_INTERNAL_OPT_SIZE": 20,
		"PRIMER_INTERNAL_MIN_SIZE": 18,
		"PRIMER_INTERNAL_MAX_SIZE": 27,
		"PRIMER_INTERNAL_OPT_TM": 60.0,  # 보통 Primer보다 Tm을 높게 잡음 (TaqMan 등)
		"PRIMER_INTERNAL_MIN_TM": 57.0,
		"PRIMER_INTERNAL_MAX_TM": 63.0,
		
		# 3. Primer 기본 조건
		"PRIMER_OPT_SIZE": 20,
		"PRIMER_MIN_SIZE": 18,
		"PRIMER_MAX_SIZE": 25,
		"PRIMER_OPT_TM": 55.0,
		"PRIMER_MIN_TM": 50.0,
		"PRIMER_MAX_TM": 60.0,
		
		# 4. 반환 개수 (Config 값 사용하거나 강제 지정)
		"PRIMER_NUM_RETURN": 3 
	}

	primer_config = app_cfg.pcr_params.primer_kwargs

	# -------------------------------------------------------------------------
	# [Step 3] Primer Design Execution
	# -------------------------------------------------------------------------
	# Probe가 들어가려면 Sequence가 어느 정도 길어야 합니다. (Amplicon 공간 확보)
	# 아래는 약 300bp 정도의 가상 서열입니다.
	seq = (
		"GATGCACCCTCCTTCTCATCCCCACCACGCCAGGCTCCTGAGCCATCCCTCTTCTTCCAGGATCCCCCTGGAACTAGTATGGAG"
	)
		
	# Target을 중간쯤 잡습니다. (Probe는 이 Target 주변이나 위에 얹혀짐)
	target_start = 52
	target_len = 1  # "GATCCCCCTGGAACTAGTAT" 부근
	use_probe = 1  # 1: Probe 포함, 0: 일반 Primer만

	if use_probe == 1:
		designer = ProbePrimerDesigner(
			template_sequence=seq, 
			target_start_index=target_start,	 
			target_end_index=target_start + target_len,		 
			config=primer_config,	   
			overrides=overrides	   
		)
	else:
		designer = PrimerDesigner(
			template_sequence=seq, 
			target_start_index=target_start,	 
			target_end_index=target_start + target_len,		 
			config=primer_config,	   
			overrides=overrides	   
		)

	try:
		results = designer.design()
		print(f"   -> Designed {len(results)} amplicons.")
	except Exception as e:
		print(f"❌ Design Failed: {e}")
		import traceback
		traceback.print_exc()
		sys.exit(1)

	if not results:
		print("⚠️ No primers found.")
		sys.exit(0)

	# -------------------------------------------------------------------------
	# [Step 4] QC 실행 (Thermo + BLAST + isPCR)
	# -------------------------------------------------------------------------
	print("\n🔍 [Step 3] Running QC Executor...")
	executor = QCExecutor(app_cfg.qc_params)
	checked_amplicons = executor.run_qc(results)

	# -------------------------------------------------------------------------
	# [Step 5] 검증 리포트
	# -------------------------------------------------------------------------
	print("\n📊 [Step 4] Validation Report")
		
	# 요약표
	df = executor.summarize(checked_amplicons)
	print(df.to_string(index=False))

	print("\n🔎 Detailed Verification:")
	for i, amp in enumerate(checked_amplicons):
		status = "✅ PASS" if amp.is_qc_pass else "❌ FAIL"
		print(f"\n[Amplicon #{i+1}] {status} | Size: {amp.product_size}bp")
		
		# 1. Probe & Primer Tm 비교 (로직 검증 핵심)
		p_tm = amp.probe.tm if amp.probe else 0.0
		f_tm = amp.forward_primer.tm
		r_tm = amp.reverse_primer.tm
		
		print(f"   [Tm Check]")
		print(f"	 🔹 Probe Tm: {p_tm:.2f} C")
		print(f"	 🔸 Fwd Tm:   {f_tm:.2f} C (Diff: {p_tm - f_tm:.2f})")
		print(f"	 🔸 Rev Tm:   {r_tm:.2f} C (Diff: {p_tm - r_tm:.2f})")
		
		# Diff가 5~10 사이인지 확인
		avg_diff = p_tm - (f_tm + r_tm)/2
		if 5.0 <= avg_diff <= 10.0 + 1.0: # 약간의 오차 허용
			print(f"	 ✅ Tm Logic Verified (Probe is ~{avg_diff:.1f}C higher)")
		else:
			print(f"	 ⚠️ Tm Logic Warning (Avg Diff: {avg_diff:.1f}C)")

		# 2. Sequence
		print(f"   [Sequences]")
		print(f"	 P: {amp.probe.sequence}")
		print(f"	 F: {amp.forward_primer.sequence}")
		print(f"	 R: {amp.reverse_primer.sequence}")

		# 3. Thermo QC 상세 (Probe 관련 QC가 돌았는지 확인)
		t_data = amp.qc_status.thermo.data
		if 'probe_hairpin_dg' in t_data:
			 print(f"   [Thermo QC]")
			 print(f"	 - Probe Hairpin dG: {t_data['probe_hairpin_dg']:.2f}")
			 print(f"	 - Hetero (F-P): {'OK' if t_data.get('hetero_fp_pass') else 'Fail'}")
			 print(f"	 - Hetero (R-P): {'OK' if t_data.get('hetero_rp_pass') else 'Fail'}")

		# 4. BLAST/isPCR 결과
		spec = amp.qc_status.specificity
		if spec.get('skipped'):
			print("   [Specificity] Skipped (Thermo Failed)")
		else:
			print(f"   [Specificity] {'OK' if spec['passed'] else 'FAIL'} (Off-targets: {spec.get('off_target_count')})")

if __name__ == "__main__":
	main()