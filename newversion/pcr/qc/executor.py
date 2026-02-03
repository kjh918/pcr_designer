import pandas as pd
from typing import List, Dict, Any, Optional

# 1. Config & Schema Import
from ..config.schema.root import AppConfig
# ✅ 핵심: 여기서 Pydantic 모델을 가져옵니다. (재정의 X)
from ..config.schema.qc import AmpliconQCStatus 

# 2. Components
from ..components.amplicon import Amplicon

# 3. QC Modules
from .thermo import ThermoChecker
from .blast import BlastSpecificityChecker
from .ispcr import IsPcrChecker

class QCExecutor:
	"""
	QC 파이프라인의 총괄 실행자.
	[최적화 전략]
	1. Thermo Check (Fast): 물리적 성질(Hairpin, Dimer) 미달 시 즉시 탈락
	2. Specificity Check (Slow): Thermo 통과한 후보만 BLAST 수행
	3. isPCR Check (Fast/Local): 최종 후보에 대해 실제 증폭 여부 검증
	"""
	def __init__(self, config: AppConfig):
		self.config = config
		self.criteria = config.qc_criteria
		
		# 각 검사기 초기화
		self.thermo_checker = ThermoChecker(self.criteria)
		self.spec_checker = BlastSpecificityChecker(config)
		
		# isPCR 검사기는 옵션에 따라 활성화
		self.ispcr_checker = None
		if self.criteria.use_ispcr_check:
			self.ispcr_checker = IsPcrChecker(config)

	def run_qc(self, amplicons: List[Amplicon]) -> List[Amplicon]:
		"""
		[Main Pipeline] 모든 앰플리콘에 대해 단계별 QC 수행
		"""
		print(f"🚀 [QC] Starting validation for {len(amplicons)} candidates...")
		
		skipped_blast = 0
		skipped_ispcr = 0

		for amp in amplicons:
			# -------------------------------------------------------------
			# Step 1: Thermo QC (매우 빠름)
			# -------------------------------------------------------------
			print(amp)
			t_res = self.thermo_checker.check(amp)
			s_res = None
			i_res = None
			final_pass = False
			# -------------------------------------------------------------
			# Step 2: Specificity QC (BLAST - 느림) -> 조건부 실행
			# -------------------------------------------------------------
			if t_res["passed"]:
				# 물성을 통과한 우량 후보만 BLAST 진행
				s_res = self.spec_checker.check_amplicon(amp)
				# ---------------------------------------------------------
				# Step 3: isPCR QC (Local Fetch - 빠름) -> 조건부 실행
				# ---------------------------------------------------------
				if s_res["passed"] and self.ispcr_checker:
					i_res = self.ispcr_checker.check_amplicon(amp)
					# 최종 통과 여부: Thermo & Blast & isPCR 모두 Pass
					final_pass = i_res["passed"]
				elif s_res["passed"] and not self.ispcr_checker:
					# isPCR 안 쓰면 BLAST 결과가 최종
					final_pass = True
					i_res = {"passed": True, "msg": "Skipped (Config)"}
				else:
					# BLAST 탈락 시 isPCR 스킵
					final_pass = False
					skipped_ispcr += 1
			else:
				# 물성 탈락 시 BLAST 스킵
				skipped_blast += 1
				s_res = {"passed": False, "skipped": True, "reason": "Thermo Failed"}
			
			# -------------------------------------------------------------
			# Step 4: 결과 저장 (Pydantic Model 사용)
			# -------------------------------------------------------------
			# 에러 원인: 여기서 Pydantic 모델(AmpliconQCStatus)을 생성해 주입해야 함
			qc_status = AmpliconQCStatus(
				thermo=t_res,
				specificity=s_res,
				ispcr=i_res,
				overall_passed=final_pass
			)
			
			# Amplicon 객체 업데이트
			amp.qc_status = qc_status	  # 객체 저장
			amp.qc_details = qc_status.model_dump() # 딕셔너리 저장 (호환성용)
			amp.is_qc_pass = final_pass
			
			# 요약 메시지 생성 (리포트용)
			amp.qc_details["fail_reason"] = self._generate_fail_msg(qc_status)

		# 로그 출력
		if skipped_blast > 0:
			print(f"   └─ Optimization: BLAST skipped for {skipped_blast} candidates (Thermo Failed).")
		
		return amplicons

	def _generate_fail_msg(self, status: AmpliconQCStatus) -> str:
		"""실패 원인을 한 줄 요약"""
		if status.overall_passed:
			return "Pass"
		
		reasons = []
		# Thermo 실패
		if status.thermo and not status.thermo.get("passed"):
			reasons.append(f"Thermo({status.thermo.get('reason', 'Fail')})")
		
		# Specificity 실패
		if status.specificity:
			if status.specificity.get("skipped"):
				pass # 이미 Thermo에서 잡힘
			elif not status.specificity.get("passed"):
				off_cnt = status.specificity.get("off_target_count", 0)
				target_ok = status.specificity.get("intended_target_found", False)
				if not target_ok:
					reasons.append("Spec(Target Not Found)")
				if off_cnt > 0:
					reasons.append(f"Spec({off_cnt} Off-targets)")

		# isPCR 실패
		if status.ispcr and not status.ispcr.get("passed"):
			 reasons.append("isPCR(Verification Failed)")

		return ", ".join(reasons) if reasons else "Unknown Fail"

	def summarize(self, amplicons: List[Amplicon]) -> pd.DataFrame:
		"""결과를 DataFrame으로 변환 (엑셀 저장용)"""
		data = []
		for amp in amplicons:
			row = {
				"ID": amp.id,
				"Status": "Pass" if amp.is_qc_pass else "Fail",
				"Reason": amp.qc_details.get("fail_reason", ""),
				
				# Basic Info
				"Chrom": amp.reference_id,
				"Size": amp.product_size,
				"Fwd_Seq": amp.forward.sequence,
				"Rev_Seq": amp.reverse.sequence,
				"Fwd_Tm": round(amp.forward.tm, 1),
				"Rev_Tm": round(amp.reverse.tm, 1),
			}
			
			# Thermo Info (Hairpin/Dimer)
			if amp.qc_status and amp.qc_status.thermo:
				t = amp.qc_status.thermo.get("data", {})
				row["Fwd_Hairpin"] = t.get("fwd_hairpin_dg", 0)
				row["Rev_Hairpin"] = t.get("rev_hairpin_dg", 0)
				row["Hetero_Dimer"] = t.get("hetero_fr_dg", 0)

			# Specificity Info
			if amp.qc_status and amp.qc_status.specificity:
				s = amp.qc_status.specificity
				row["Off_Targets"] = s.get("off_target_count", 0)
				
			data.append(row)
			
		return pd.DataFrame(data)