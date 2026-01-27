from typing import List, Dict, Any
import pandas as pd

from ..config.schema.qc import QCParams
from ..components import Amplicon

# Checkers
from .thermo import ThermoChecker
from .blast import BlastSpecificityChecker
from .types import QCResult

class AmpliconQCStatus:
	"""Amplicon의 QC 상세 결과를 담는 객체"""
	def __init__(self, is_pass: bool, thermo: QCResult, specificity: Dict[str, Any]):
		self.is_pass = is_pass
		self.thermo = thermo
		self.specificity = specificity # BLAST + isPCR 결과 Dict

	def to_summary(self) -> Dict[str, Any]:
		"""Summary Report용 요약 데이터 변환"""
		# Thermo 실패로 인해 BLAST가 스킵되었는지 확인
		is_spec_skipped = self.specificity.get("skipped", False)
		
		spec_status = "Skipped" if is_spec_skipped else ("OK" if self.specificity['passed'] else "Fail")
		
		off_targets = "N/A"
		if not is_spec_skipped:
			off_targets = self.specificity['off_target_count']

		target_found = "-"
		if not is_spec_skipped:
			target_found = "Yes" if self.specificity['intended_target_found'] else "No"

		return {
			"Final": "PASS" if self.is_pass else "FAIL",
			"Thermo": "OK" if self.thermo.passed else "Fail",
			"Specificity": spec_status,
			"Off_Targets": off_targets,
			"Target_Found": target_found
		}

class QCExecutor:
	def __init__(self, qc_params: QCParams):
		self.params = qc_params
		
		# 1. Thermo Checker (물리적 성질)
		self.thermo_checker = ThermoChecker(qc_params)
		
		# 2. Specificity Checker (BLAST Search + isPCR Verify)
		self.spec_checker = BlastSpecificityChecker(qc_params)

	def run_qc(self, amplicons: List[Amplicon]) -> List[Amplicon]:
		"""
		[Main Pipeline - Optimized]
		Thermo 통과한 경우에만 BLAST/isPCR 수행 (속도 최적화)
		"""
		print(f"🚀 [QC Executor] Starting pipeline for {len(amplicons)} amplicons...")
		
		skipped_count = 0

		for amp in amplicons:
			# ----------------------------------------------------------------
			# 1. Thermo QC 실행 (빠름)
			# ----------------------------------------------------------------
			t_res = self.thermo_checker.check(amp) # returns QCResult

			# ----------------------------------------------------------------
			# 2. Specificity QC 실행 (BLAST: 느림) -> 조건부 실행
			# ----------------------------------------------------------------
			if t_res.passed:
				# Thermo를 통과했으므로 비싼 BLAST 연산 수행
				s_res = self.spec_checker.check_amplicon(amp)
				final_pass = s_res['passed']
			else:
				# Thermo 탈락 -> BLAST 스킵 (Dummy Result 생성)
				skipped_count += 1
				s_res = {
					"passed": False,
					"skipped": True, # 스킵 마커
					"off_target_count": -1,
					"intended_target_found": False,
					"off_targets": [],
					"alignment": None
				}
				final_pass = False

			# ----------------------------------------------------------------
			# 3. 결과 객체 생성 및 주입
			# ----------------------------------------------------------------
			status = AmpliconQCStatus(
				is_pass=final_pass,
				thermo=t_res,
				specificity=s_res
			)
			
			# Amplicon 객체에 속성 주입
			amp.qc_status = status	   
			amp.is_qc_pass = final_pass
			
			# 편의 속성 (스킵된 경우 0이나 -1 처리)
			amp.off_target_count = s_res.get('off_target_count', 0)

		print(f"   -> Optimization: Skipped BLAST for {skipped_count} amplicons (Thermo Failed).")
		return amplicons

	def summarize(self, amplicons: List[Amplicon]) -> pd.DataFrame:
		"""QC 결과 요약 DataFrame 생성"""
		data = []
		for i, amp in enumerate(amplicons):
			if not hasattr(amp, 'qc_status'): continue
			
			row = {
				"ID": f"Amp_{i+1}",
				"Size": amp.product_size,
				"Tm": round(amp.tm, 2)
			}
			# 상세 결과 병합
			row.update(amp.qc_status.to_summary())
			
			# 실패 원인 요약
			fails = []
			if not amp.qc_status.thermo.passed: 
				fails.append("Thermo")
			elif not amp.qc_status.specificity.get('passed'): 
				# Thermo는 통과했는데 Spec에서 떨어진 경우만 기록
				fails.append("Specificity")
			
			row["Fail_Reason"] = ", ".join(fails) if fails else "-"
			data.append(row)
			
		return pd.DataFrame(data)