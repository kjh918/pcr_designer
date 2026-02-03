import sys
import os
import pandas as pd

# 프로젝트 루트 경로 추가 (모듈 import를 위해)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pcr.config.loader import ConfigLoader
from pcr.config.schema.root import BaseDesignInput
from pcr.designers.qpcr import QPCRPrimerDesigner  # 작성하신 클래스

def run_test():
	print("🧪 [Test] Starting qPCR Primer Design Test...\n")

	# ---------------------------------------------------------
	# 1. 설정 로드
	# ---------------------------------------------------------
	try:
		loader = ConfigLoader()
		# default.yaml에 probe_kwargs 등이 정의되어 있어야 함
		config = loader.load(preset_name="default")
		print("✅ Config loaded.")
		
		# 확인: 설정된 Tm Diff 출력
		criteria = config.qc_criteria.probe
		print(f"   ℹ️  Target Tm Diff: Min {criteria.min_primer_probe_tm_diff}°C ~ Max {criteria.max_primer_probe_tm_diff}°C")
		
	except Exception as e:
		print(f"❌ Config Load Error: {e}")
		return

	# ---------------------------------------------------------
	# 2. 테스트 데이터 준비 (EGFR Exon 21 L858R 주변 예시)
	# ---------------------------------------------------------
	# 약 400bp 길이의 DNA 서열
	dummy_seq = (
		"AGCCTCTTACACCCAGTGGAGAAGCTCCCAACCAAGCCAACAGGTCCTGGAGCTTTGGGGCCACTCTACAAACAAAAAAC"
		"AAAAACAAAAACAAAAACAAAAACAAAAAACCAAGCCAACAGGTCCTGGAGCTTTGGGGCCACTCTACAAACAAACAAAA"
		"GCAGGGGTTGGCCCTGCCCAACCAAGCCAACTGGTCCTGGAGCTTTGGGGCCACTCTACAAACAAACAAAACCAAGCCAA"
		"CAGGTCCTGGAGCTTTGGGGCCACTCTACAAACAAACAAAACCAAGCCAACAGGTCCTGGAGCTTTGGGGCCACTCTACA"
		"AACAAACAAAACCAAGCCAACAGGTCCTGGAGCTTTGGGGCCACTCTACAAACAAACAAAA"
	)
	# 실제로는 Primer3가 잘 찾을 수 있는 조금 더 랜덤하고 복잡한 서열이 좋지만, 테스트용으로 사용
	# 타겟을 중앙에 잡음
	target_start = 150
	target_end = 151 # SNP 위치라고 가정 (1bp)

	print(f"   ℹ️  Sequence Length: {len(dummy_seq)}bp")
	print(f"   ℹ️  Target Region: {target_start}-{target_end}")

	# ---------------------------------------------------------
	# 3. Input 객체 생성
	# ---------------------------------------------------------
	design_input = BaseDesignInput(
		name="Test_L858R",
		template_sequence=dummy_seq,
		target_start=target_start,
		target_end=target_end,
		config=config,
		reference_name="Synthesized",
		overrides={
			"PRIMER_NUM_RETURN": 100  # 테스트니까 Probe 3개까지만 진행
		}
	)

	# ---------------------------------------------------------
	# 4. Designer 실행 (QPCR 모드)
	# ---------------------------------------------------------
	print("\n🚀 Running QPCRPrimerDesigner...")
	designer = QPCRPrimerDesigner(design_input)
		
	# design() 메서드 호출 (QC 없음, 순수 디자인)
	output = designer.design()
	print(output)
	# ---------------------------------------------------------
	# 5. 결과 검증 및 출력
	# ---------------------------------------------------------
	if output.status != "success":
		print(f"❌ Design Failed: {output.status}")
		if output.error_msg:
			print(f"   Error: {output.error_msg}")
		if output.log_messages:
			print(f"   Logs: {output.log_messages}")
		return

	amplicons = output.amplicons
	print(f"\n✅ Design Success! Found {len(amplicons)} candidates.\n")

	# 결과 상세 출력 (데이터프레임 형태)
	data = []
	for amp in amplicons:
		# Tm 차이 검증
		tm_diff_fwd = amp.probe.tm - amp.forward.tm
		tm_diff_rev = amp.probe.tm - amp.reverse.tm
		
		row = {
			"ID": amp.id,
			"Probe_Tm": round(amp.probe.tm, 2),
			"Fwd_Tm": round(amp.forward.tm, 2),
			"Rev_Tm": round(amp.reverse.tm, 2),
			"Diff_F (P-F)": round(tm_diff_fwd, 2), # 양수여야 함 (Probe > Primer)
			"Diff_R (P-R)": round(tm_diff_rev, 2), # 양수여야 함
			"Product_Size": amp.product_size,
			"Probe_Seq": amp.probe.sequence
		}
		data.append(row)

	df = pd.DataFrame(data)
		
	# 터미널에 예쁘게 출력
	pd.set_option('display.max_columns', None)
	pd.set_option('display.width', 1000)
	print(df.to_string(index=False))

	print("\n---------------------------------------------------------")
	print("🧪 검증 포인트:")
	print("1. ID가 'Test_L858R_Px_y' 형식인가? (Probe-First 구조)")
	print("2. Diff_F와 Diff_R이 설정범위(예: 5~10도) 내에 있는가?")
	print("3. Probe_Seq가 G로 시작하지 않는가? (필터링 확인)")
	print("---------------------------------------------------------")

if __name__ == "__main__":
	run_test()