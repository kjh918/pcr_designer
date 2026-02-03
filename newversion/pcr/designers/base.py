"""
pcr/designers/base.py
PCR Primer Design의 기본 클래스.
(QC 및 Ranking 로직 제거 -> Primer3 원본 결과 반환)
"""
import primer3
from typing import Dict, Any, List

# 1. Config & Schema
from ..config.schema.root import BaseDesignInput, BaseDesignOutput
from ..components.primer import Primer, Probe
from ..components.amplicon import Amplicon

# Ranker, QCExecutor import 모두 제거됨

class BasePrimerDesigner:
	"""
	Standard PCR Designer (Raw Output Version)
	Input -> Primer3 -> Object Conversion -> Output
	"""
	ASSAY_TYPE = "Base-PCR"

	def __init__(self, input_data: BaseDesignInput):
		self.input = input_data
		self.config = input_data.config
		
		# Assay 이름 설정
		self.design_name = f"{self.input.name}_{self.ASSAY_TYPE}"
		
		# Ranker, QCExecutor 초기화 로직 삭제됨

		# Primer3 설정 컨테이너
		self.seq_args: Dict[str, Any] = {}
		self.global_args: Dict[str, Any] = {}
		
		self._prepare_primer3_args()

	def _prepare_primer3_args(self):
		"""Primer3 입력 인자 준비"""
		target_len = self.input.target_end - self.input.target_start
		
		self.seq_args = {
			'SEQUENCE_ID': self.design_name,
			'SEQUENCE_TEMPLATE': self.input.template_sequence,
			'SEQUENCE_TARGET': [self.input.target_start, target_len]
		}

		# Config -> Primer3 Global Args 변환
		self.global_args = self.config.pcr_params.primer_kwargs.to_global_args()

		# Probe 설정 (옵션)
		if self.config.pcr_params.probe_kwargs:
			probe_args = self.config.pcr_params.probe_kwargs.to_global_args()
			self.global_args.update(probe_args)
			self.global_args['PRIMER_PICK_INTERNAL_OLIGO'] = 1
		else:
			self.global_args['PRIMER_PICK_INTERNAL_OLIGO'] = 0

		# Runtime Overrides
		if self.input.overrides:
			self.global_args.update(self.input.overrides)

	def design(self) -> BaseDesignOutput:
		"""
		[Main Pipeline]
		Primer3 실행 -> 결과 객체 변환 -> 반환
		"""
		try:
			# Step 1: Primer3 실행
			raw_result = primer3.bindings.design_primers(self.seq_args, self.global_args)
			
			# Step 2: 객체 변환
			all_candidates = self._process_results(raw_result)
			
			if not all_candidates:
				return BaseDesignOutput(
					amplicons=[], 
					status="no_candidates",
					log_messages=["Primer3 found 0 candidates."]
				)

			# Step 3: 결과 반환
			# 별도의 Ranking 없이, Primer3가 준 순서(Penalty 낮은 순) 그대로 반환
			# 단, 너무 많을 수 있으므로 설정된 top_k 만큼 자르긴 함 (원치 않으면 이 부분도 제거 가능)
			limit = self.input.overrides.get("return_top_k", 500) 
			final_list = all_candidates[:limit]
			
			return BaseDesignOutput(
				amplicons=final_list,
				status="success",
				log_messages=[
					f"Primer3 Returned: {len(all_candidates)}",
					f"Returned: {len(final_list)}",
					"Pipeline: Design Only (No QC/Ranking)"
				]
			)

		except Exception as e:
			import traceback
			traceback.print_exc()
			return BaseDesignOutput(
				amplicons=[], 
				status="error", 
				error_msg=str(e)
			)

	def _process_results(self, result: Dict[str, Any]) -> List[Amplicon]:
		"""Primer3 결과를 Amplicon 객체 리스트로 변환"""
		num_returned = result.get('PRIMER_PAIR_NUM_RETURNED', 0)
		amplicons = []
		
		for i in range(num_returned):
			fwd = Primer.from_primer3(result, i, "LEFT")
			rev = Primer.from_primer3(result, i, "RIGHT")
			
			if not fwd or not rev: continue
			probe = None
			if self.config.pcr_params.probe_kwargs and f"PRIMER_INTERNAL_{i}_SEQUENCE" in result:
				probe = Probe.from_primer3(result, i)

			pair_penalty = float(result.get(f'PRIMER_PAIR_{i}_PENALTY', 0.0))
			
			amp_id = f"{self.input.name}_{i}"
			
			amp = Amplicon(
				id=amp_id,
				forward=fwd,
				reverse=rev,
				probe=probe,
				template_sequence=self.input.template_sequence,
				target_start_index=self.input.target_start,
				target_end_index=self.input.target_end,
				reference_id=self.input.reference_name,
				pair_penalty=pair_penalty
			)
			
			# QC 패스 여부는 나중에 외부에서 판단하므로 일단 True(또는 False)로 둠
			amp.is_qc_pass = True 
			
			amplicons.append(amp)
			
		return amplicons