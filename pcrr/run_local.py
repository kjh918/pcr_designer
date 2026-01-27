import sys
import os
import pandas as pd

# 1. src 폴더를 라이브러리 경로에 추가
sys.path.append(os.path.abspath("src"))

from src.pcr.config.loader import load_config
from src.pcr.designers.base import BasePrimerDesigner
# ✅ Executor 임포트 (개별 Checker 임포트 불필요)
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
	# [Step 2] Primer Design (Simulation)
	# -------------------------------------------------------------------------
	# 테스트용 짧은 서열 (실제 구동 시에는 유효한 유전자 서열 사용 권장)
	seq = "GATGCACCCTCCTTCTCATCCCCACCACGCCAGGCTCCTGAGCCATCCCTCTTCTTCCAGGATCCCCCTGGAACTAGTATGGAG" 
		
	# 일부러 GC 조건을 넓혀서 디자인이 되도록 유도
	overrides = {
		"PRIMER_OPT_GC_PERCENT": 50.0,
		"PRIMER_MIN_GC": 20.0,
		"PRIMER_MAX_GC": 80.0,
		"PRIMER_NUM_RETURN": 10
	}

	primer_config = app_cfg.pcr_params.primer_kwargs
	designer = BasePrimerDesigner(
		template_sequence=seq, 
		target_start_index=51,	 
		target_end_index=51,		 
		config=primer_config,	   
		overrides=overrides	   
	)
	print("\n🧪 [Step 2] Designing Primers...")
	try:
		results = designer.design()
		print(f"   -> Designed {len(results)} amplicons.")
	except Exception as e:
		print(f"❌ Design Failed: {e}")
		sys.exit(1)

	if not results:
		print("⚠️ No primers found. Exiting.")
		sys.exit(0)

	# -------------------------------------------------------------------------
	# [Step 3] QC Execution (Updated)
	# -------------------------------------------------------------------------
	print("\n🔍 [Step 3] Running QC Executor (Thermo + BLAST + isPCR)...")
		
	# Executor 인스턴스 생성 및 실행
	executor = QCExecutor(app_cfg.qc_params)
	checked_amplicons = executor.run_qc(results)

	# -------------------------------------------------------------------------
	# [Step 4] 결과 리포팅
	# -------------------------------------------------------------------------
	print("\n📊 [Step 4] QC Summary Report")
		
	# 1. DataFrame 요약 출력
	df = executor.summarize(checked_amplicons)
	print(df.to_string(index=False))

	# 2. 상세 정보 출력 (Pass된 것 위주)
	print("\n🔎 Detailed View:")
	for i, amp in enumerate(checked_amplicons):
		status_icon = "✅ PASS" if amp.is_qc_pass else "❌ FAIL"
		print(f"\n[Amplicon #{i+1}] {status_icon}")
		
		# Specificity 결과 접근
		spec_res = amp.qc_status.specificity
		
		if spec_res['passed']:
			print(f"   - Specificity: OK (Unique Target Confirmed)")
		else:
			print(f"   - Specificity: FAIL (Found {spec_res['off_target_count']} off-targets)")

		# Alignment 시각화 (Target 찾았을 경우)
		if spec_res.get('alignment'):
			print(f"   - Target Alignment ({spec_res['alignment'].get('chrom', 'N/A')}):")
			viz = spec_res['alignment'].get('visualization', '')
			for line in viz.split('\n'):
				print(f"	 {line}")
		elif not spec_res.get('intended_target_found'):
			 print("   - ⚠️ Intended target not found in Genome (Check Sequence/DB).")

if __name__ == "__main__":
	main()