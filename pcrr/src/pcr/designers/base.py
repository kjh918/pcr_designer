import primer3
from typing import Any, Dict, List, Optional
from pydantic import BaseModel, Field

# 프로젝트 구조에 맞는 Import
from ..components import Amplicon, Primer
from ..config.schema.pcr import PrimerKwargs

# ----------------------------------------------------------------------
# [Base] 공통 부모 클래스 (기존 구조 유지)
# ----------------------------------------------------------------------
class BasePrimerDesigner:
	"""
	Primer3 실행을 위한 공통 기반 클래스
	Config 객체와 Overrides 딕셔너리를 받아 설정을 초기화함.
	"""
	def __init__(
		self,
		template_sequence: str,
		target_start_index: int,
		target_end_index: int,
		config: PrimerKwargs,
		overrides: Optional[Dict[str, Any]] = None,
	) -> None:
		self.template_sequence = template_sequence
		self.target_start_index = int(target_start_index)
		self.target_end_index = int(target_end_index)
		
		# ✅ 기존 구조대로 Config와 Overrides 저장
		self.cfg = config
		self.overrides = overrides or {}

		# Primer3용 Args 초기화
		self.primer3_seq_args: Dict[str, Any] = {}
		self.primer3_global_args: Dict[str, Any] = {}
		
		self.primer3_result: Optional[Dict[str, Any]] = None
		self.amplicon_list: List[Amplicon] = []

		# 초기 설정 로드
		self._init_primer3_args()

	def _init_primer3_args(self) -> None:
		"""기본 Sequence 및 Global Args 설정"""
		target_len = self.target_end_index - self.target_start_index + 1
		
		self.primer3_seq_args = {
			"SEQUENCE_ID": "generic_design",
			"SEQUENCE_TEMPLATE": self.template_sequence,
			"SEQUENCE_TARGET": [self.target_start_index, target_len],
		}

		# Config -> Primer3 Args 변환
		if hasattr(self.cfg, "to_global_args"):
			self.primer3_global_args = self.cfg.to_global_args()
		else:
			self.primer3_global_args = self.cfg.to_primer3_args()

		# Overrides 적용 (Config보다 우선순위 높음)
		if self.overrides:
			self.primer3_global_args.update(self.overrides)

	def run_primer3(self) -> None:
		"""설정된 Args로 Primer3 실행"""
		try:
			# 최신 함수명 사용
			self.primer3_result = primer3.bindings.design_primers(
				seq_args=self.primer3_seq_args,
				global_args=self.primer3_global_args,
			)
		except Exception as e:
			# 디자인 실패 시 빈 결과 처리 혹은 에러 로깅
			self.primer3_result = {}
			# raise RuntimeError(f"Primer3 Execution Failed: {e}") 

	def design(self) -> List[Amplicon]:
		"""[Public API] 디자인 실행"""
		self.run_primer3()
		return self._build_amplicons()

	def _build_amplicons(self) -> List[Amplicon]:
		raise NotImplementedError("Subclasses must implement _build_amplicons")
		
	def reset(self) -> None:
		"""설정 초기화 (재사용 시)"""
		self._init_primer3_args()
		self.primer3_result = None
		self.amplicon_list = []


# ----------------------------------------------------------------------
# [PrimerDesigner] 일반 프라이머 (Probe 없음)
# ----------------------------------------------------------------------
class PrimerDesigner(BasePrimerDesigner):
	"""일반 Primer Pair 디자인 클래스"""
		
	def _build_amplicons(self) -> List[Amplicon]:
		if not self.primer3_result:
			return []

		n_pairs = int(self.primer3_result.get("PRIMER_PAIR_NUM_RETURNED", 0))
		amplicons: List[Amplicon] = []

		for i in range(n_pairs):
			# Forward
			fwd = Primer(
				template_sequence=self.template_sequence,
				sequence=self.primer3_result.get(f"PRIMER_LEFT_{i}_SEQUENCE"),
				tm=self.primer3_result.get(f"PRIMER_LEFT_{i}_TM"),
				strand="forward",
				primer_type="forward_primer",
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				binding_start_index=self.primer3_result.get(f"PRIMER_LEFT_{i}")[0],
				binding_end_index=self.primer3_result.get(f"PRIMER_LEFT_{i}")[0] + self.primer3_result.get(f"PRIMER_LEFT_{i}")[1]
			)
			
			# Reverse
			rev_data = self.primer3_result.get(f"PRIMER_RIGHT_{i}")
			rev_start, rev_len = rev_data[0], rev_data[1]
			rev = Primer(
				template_sequence=self.template_sequence,
				sequence=self.primer3_result.get(f"PRIMER_RIGHT_{i}_SEQUENCE"),
				tm=self.primer3_result.get(f"PRIMER_RIGHT_{i}_TM"),
				strand="reverse",
				primer_type="reverse_primer",
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				binding_start_index=rev_start - rev_len + 1,
				binding_end_index=rev_start + 1
			)

			amp = Amplicon(
				template_sequence=self.template_sequence,
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				forward_primer=fwd,
				reverse_primer=rev,
				product_size=self.primer3_result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
				tm=self.primer3_result.get(f"PRIMER_PAIR_{i}_PRODUCT_TM")
			)
			amplicons.append(amp)

		return amplicons


# ----------------------------------------------------------------------
# [ProbePrimerDesigner] Probe 우선 디자인 -> Tm 조정 -> Primer 디자인
# ----------------------------------------------------------------------
class ProbePrimerDesigner(BasePrimerDesigner):
	"""
	[Advanced Logic]
	1. Overrides 설정을 읽어 Probe만 먼저 디자인 (Primer Off)
	2. 디자인된 Probe들의 실제 Tm을 확인
	3. 각 Probe Tm에 맞춰 Primer Tm Window를 동적으로 설정하여 Pair 디자인
	"""
	# ✅ __init__은 BasePrimerDesigner 것을 그대로 사용하므로 재정의 불필요

	def _configure_probe_only(self) -> None:
		"""Probe만 뽑도록 설정 변경"""
		
		self.primer3_global_args["PRIMER_PICK_LEFT_PRIMER"] = 0
		self.primer3_global_args["PRIMER_PICK_RIGHT_PRIMER"] = 0
		self.primer3_global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 1
		# 3. Probe 개수 및 조건 (Overrides에서 가져오기)
		# run_local.py의 overrides가 이미 update되어 있으므로, 필요한 것만 확인
		n_return = self.overrides.get("PRIMER_NUM_RETURN", 5)
		self.primer3_global_args["PRIMER_INTERNAL_NUM_RETURN"] = n_return
		#self.primer3_global_args["PRIMER_NUM_RETURN"] = 0 

	def _build_probe_only_amplicons(self) -> List[Amplicon]:
		"""Probe만 있는 임시 Amplicon 리스트 생성"""
		if not self.primer3_result:
			return []
			
		n_probes = int(self.primer3_result.get("PRIMER_INTERNAL_NUM_RETURNED", 0))
		amplicons = []

		for i in range(n_probes):
			seq = self.primer3_result.get(f"PRIMER_INTERNAL_{i}_SEQUENCE")
			tm_val = self.primer3_result.get(f"PRIMER_INTERNAL_{i}_TM")
			
			if not seq: continue

			probe = Primer(
				template_sequence=self.template_sequence,
				sequence=seq,
				tm=float(tm_val) if tm_val else 0.0,
				strand="forward",
				primer_type="probe",
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
			)
			
			# Primer 없는 Amplicon 생성
			amplicons.append(Amplicon(
				template_sequence=self.template_sequence,
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				forward_primer=None,
				reverse_primer=None,
				probe=probe
			))
		return amplicons

	def _design_primers_for_probe(self, probe_amp: Amplicon) -> List[Amplicon]:
		"""고정된 Probe에 맞춰 Primer Pair 디자인"""
		probe = probe_amp.probe
		if not probe or not probe.tm:
			return []

		# 1. 설정 초기화 (Primer On)
		self.reset()
		self.primer3_global_args["PRIMER_PICK_LEFT_PRIMER"] = 1
		self.primer3_global_args["PRIMER_PICK_RIGHT_PRIMER"] = 1
		self.primer3_global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 1 
		
		# 2. Probe 고정 (Primer3에게 이 Probe를 쓰라고 강제)
		self.primer3_seq_args["SEQUENCE_INTERNAL_OLIGO"] = probe.sequence
		
		# 3. ✅ Tm Window 동적 설정 (Probe Tm - Diff)
		# overrides에서 설정값 가져오기 (없으면 기본값)
		min_diff = self.overrides.get("min_primer_probe_tm_diff", 5.0)
		max_diff = self.overrides.get("max_primer_probe_tm_diff", 10.0)
		
		target_tm = probe.tm
		
		self.primer3_global_args["PRIMER_OPT_TM"] = target_tm - min_diff
		self.primer3_global_args["PRIMER_MAX_TM"] = target_tm - min_diff + 2.0 
		self.primer3_global_args["PRIMER_MIN_TM"] = target_tm - max_diff
		
		# 4. 실행
		self.run_primer3()
		
		# 5. 결과 파싱
		if not self.primer3_result:
			return []

		n_pairs = int(self.primer3_result.get("PRIMER_PAIR_NUM_RETURNED", 0))
		results = []
		
		for i in range(n_pairs):
			# Fwd
			fwd = Primer(
				template_sequence=self.template_sequence,
				sequence=self.primer3_result.get(f"PRIMER_LEFT_{i}_SEQUENCE"),
				tm=self.primer3_result.get(f"PRIMER_LEFT_{i}_TM"),
				strand="forward", primer_type="forward_primer",
				target_start_index=self.target_start_index, target_end_index=self.target_end_index,
				binding_start_index=self.primer3_result.get(f"PRIMER_LEFT_{i}")[0],
				binding_end_index=self.primer3_result.get(f"PRIMER_LEFT_{i}")[0] + self.primer3_result.get(f"PRIMER_LEFT_{i}")[1]
			)
			# Rev
			rev_data = self.primer3_result.get(f"PRIMER_RIGHT_{i}")
			rev = Primer(
				template_sequence=self.template_sequence,
				sequence=self.primer3_result.get(f"PRIMER_RIGHT_{i}_SEQUENCE"),
				tm=self.primer3_result.get(f"PRIMER_RIGHT_{i}_TM"),
				strand="reverse", primer_type="reverse_primer",
				target_start_index=self.target_start_index, target_end_index=self.target_end_index,
				binding_start_index=rev_data[0] - rev_data[1] + 1,
				binding_end_index=rev_data[0] + 1
			)
			
			# Probe 포함된 최종 Amplicon
			results.append(Amplicon(
				template_sequence=self.template_sequence,
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				forward_primer=fwd,
				reverse_primer=rev,
				probe=probe, # 고정된 Probe
				product_size=self.primer3_result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
				tm=self.primer3_result.get(f"PRIMER_PAIR_{i}_PRODUCT_TM")
			))
			
		return results

	def design(self) -> List[Amplicon]:
		"""[Override] Probe First Strategy 실행"""
		
		# Step 1: Probe Only Design
		self.reset()
		self._configure_probe_only()
		self.run_primer3()
		probe_candidates = self._build_probe_only_amplicons()
		
		final_results = []
		
		# Step 2: Loop through Probes & Design Primers
		for probe_amp in probe_candidates:
			paired_amplicons = self._design_primers_for_probe(probe_amp)
			final_results.extend(paired_amplicons)
			
		self.amplicon_list = final_results
		return final_results

	def _build_amplicons(self) -> List[Amplicon]:
		# design() 메서드에서 직접 로직을 제어하므로 사용하지 않음
		return []