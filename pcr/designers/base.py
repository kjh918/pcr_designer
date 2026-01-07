from __future__ import annotations

from typing import Any, Dict, List, Optional

import primer3
from pcr.components import Primer, Amplicon

# ----------------------------------------------------------------------
# Base
# ----------------------------------------------------------------------
class BasePrimerDesigner:
	"""
	primer3를 호출하기 위한 공통 베이스 클래스.
	- primer3 인자 구성
	- run_primer3
	- Amplicon 생성은 서브클래스에서 구현
	"""

	DEFAULT_SALT_MONOVALENT: float = 50.0
	DEFAULT_SALT_DIVALENT: float = 1.5
	DEFAULT_DNTP_CONC: float = 0.6
	DEFAULT_DNA_CONC: float = 50.0

	def _init_primer3_args(self) -> None:
		"""primer3 인자들을 초기 상태로 세팅."""
		self.primer3_seq_args = {
			"SEQUENCE_ID": "PRIMER",
			"SEQUENCE_TEMPLATE": self.template_sequence,
		}
		self.primer3_global_args = {
			"PRIMER_TASK": "generic",
			"PRIMER_NUM_RETURN": self.n_primers,
			"PRIMER_PICK_LEFT_PRIMER": int(self.forward_primer),
			"PRIMER_PICK_RIGHT_PRIMER": int(self.reverse_primer),
			"PRIMER_PICK_INTERNAL_OLIGO": 0,
		}
		self._configure_primer_common()

	def __init__(
		self,
		template_sequence: str,
		target_start_index: int,
		target_end_index: int,
		*,
		min_amplicon_length: int = 80,
		max_amplicon_length: int = 100,
		n_primers: int = 10,
		max_tm_difference: float = 2.0,
		forward_primer: bool = True,
		reverse_primer: bool = True,
		opt_length: int = 25,
		min_length: int = 20,
		max_length: int = 30,
		opt_tm: float = 60.0,
		min_tm: float = 50.0,
		max_tm: float = 70.0,
		opt_gc: float = 45.0,
		min_gc: float = 35.0,
		max_gc: float = 65.0,
		reference_template_sequence: Optional[str] = None,
		primer3_seq_args: Optional[Dict[str, Any]] = None,
		primer3_global_args: Optional[Dict[str, Any]] = None,
	) -> None:
		# 기본 서열 정보
		self.template_sequence: str = template_sequence
		self.reference_template_sequence: str = reference_template_sequence or template_sequence
		self.target_start_index: int = target_start_index
		self.target_end_index: int = target_end_index

		# 제품 길이/프라이머 조건
		self.min_amplicon_length: int = min_amplicon_length
		self.max_amplicon_length: int = max_amplicon_length
		self.n_primers: int = n_primers
		self.max_tm_difference: float = max_tm_difference

		self.forward_primer: bool = forward_primer
		self.reverse_primer: bool = reverse_primer

		self.opt_length: int = opt_length
		self.min_length: int = min_length
		self.max_length: int = max_length

		self.opt_tm: float = opt_tm
		self.min_tm: float = min_tm
		self.max_tm: float = max_tm

		self.opt_gc: float = opt_gc
		self.min_gc: float = min_gc
		self.max_gc: float = max_gc

		# primer3 인자
		self.primer3_seq_args: Dict[str, Any] = {
			"SEQUENCE_ID": "PRIMER",
			"SEQUENCE_TEMPLATE": self.template_sequence,
		}
		self.primer3_global_args: Dict[str, Any] = {
			"PRIMER_TASK": "generic",
			"PRIMER_NUM_RETURN": self.n_primers,
			"PRIMER_PICK_LEFT_PRIMER": int(self.forward_primer),
			"PRIMER_PICK_RIGHT_PRIMER": int(self.reverse_primer),
			"PRIMER_PICK_INTERNAL_OLIGO": 0,  # 기본은 probe 없음
		}

		# 공통 primer 조건 세팅
		self._configure_primer_common()

		# 사용자 커스텀 인자 merge
		if primer3_seq_args:
			self.update_primer3_seq_args(primer3_seq_args)
		if primer3_global_args:
			self.update_primer3_global_args(primer3_global_args)

		# 결과
		self.primer3_result: Optional[Dict[str, Any]] = None
		self.amplicon_list: List[Amplicon] = []

	# ------------------------------------------------------------------
	# 공통 설정
	# ------------------------------------------------------------------
	def _configure_primer_common(self) -> None:
		target_len = self.target_end_index - self.target_start_index + 1

		self.update_primer3_seq_args(
			{
				"SEQUENCE_TARGET": [self.target_start_index, target_len],
			}
		)

		self.update_primer3_global_args(
			{
				"PRIMER_PAIR_MAX_DIFF_TM": self.max_tm_difference,
				"PRIMER_OPT_SIZE": self.opt_length,
				"PRIMER_MIN_SIZE": self.min_length,
				"PRIMER_MAX_SIZE": self.max_length,
				"PRIMER_OPT_TM": self.opt_tm,
				"PRIMER_MIN_TM": self.min_tm,
				"PRIMER_MAX_TM": self.max_tm,
				"PRIMER_OPT_GC_PERCENT": self.opt_gc,
				"PRIMER_MIN_GC": self.min_gc,
				"PRIMER_MAX_GC": self.max_gc,
				"PRIMER_PRODUCT_SIZE_RANGE": [
					self.min_amplicon_length,
					self.max_amplicon_length,
				],
			}
		)

	# ------------------------------------------------------------------
	# primer3 인자 업데이트
	# ------------------------------------------------------------------
	def update_primer3_seq_args(self, args: Dict[str, Any]) -> None:
		self.primer3_seq_args.update(args)

	def update_primer3_global_args(self, args: Dict[str, Any]) -> None:
		self.primer3_global_args.update(args)

	# ------------------------------------------------------------------
	# primer3 실행 및 결과 → Amplicon 빌드 (템플릿 메서드)
	# ------------------------------------------------------------------
	def run_primer3(self) -> None:
		self.primer3_result = primer3.bindings.designPrimers(
			seq_args=self.primer3_seq_args,
			global_args=self.primer3_global_args,
		)
		self.amplicon_list = self._build_amplicons()

	def _build_amplicons(self) -> List[Amplicon]:
		"""
		서브클래스에서 구현:
		primer3_result를 Amplicon 리스트로 변환.
		"""
		raise NotImplementedError

	def design(self) -> List[Amplicon]:
		self.run_primer3()
		return self.amplicon_list


	def reset(self) -> None:
		"""
		한 번 run_primer3()를 돌린 뒤,
		다른 조건/타겟으로 다시 쓰고 싶을 때 내부 상태를 초기화.
		"""
		self._init_primer3_args()
		self.primer3_result = None
		self.amplicon_list = []
# ----------------------------------------------------------------------
# Primer only
# ----------------------------------------------------------------------
class PrimerDesigner(BasePrimerDesigner):
	"""
	forward / reverse primer만 디자인하는 클래스.
	"""

	def _build_amplicons(self) -> List[Amplicon]:
		assert self.primer3_result is not None

		n_forward = self.primer3_result.get("PRIMER_LEFT_NUM_RETURNED", 0)
		n_reverse = self.primer3_result.get("PRIMER_RIGHT_NUM_RETURNED", 0)
		n_pairs = self.primer3_result.get("PRIMER_PAIR_NUM_RETURNED", 0)

		n_designed = max(n_forward, n_reverse, n_pairs)
		amplicons: List[Amplicon] = []

		for rank in range(n_designed):
			forward: Optional[Primer] = None
			reverse: Optional[Primer] = None

			if self.primer3_result.get(f"PRIMER_LEFT_{rank}") is not None:
				forward = Primer(
					template_sequence=self.template_sequence,
					reference_template_sequence=self.reference_template_sequence,
					sequence=self.primer3_result[f"PRIMER_LEFT_{rank}_SEQUENCE"],
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					strand="forward",
					primer_type="forward",
				)

			if self.primer3_result.get(f"PRIMER_RIGHT_{rank}") is not None:
				reverse = Primer(
					template_sequence=self.template_sequence,
					reference_template_sequence=self.reference_template_sequence,
					sequence=self.primer3_result[f"PRIMER_RIGHT_{rank}_SEQUENCE"],
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					strand="reverse",
					primer_type="reverse",
				)

			amplicon = Amplicon(
				template_sequence=self.template_sequence,
				reference_template_sequence=self.reference_template_sequence,
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				forward_primer=forward,
				reverse_primer=reverse,
			)
			amplicons.append(amplicon)

		return amplicons


# ----------------------------------------------------------------------
# Primer + Probe
# ----------------------------------------------------------------------
class ProbePrimerDesigner(BasePrimerDesigner):
	"""
	primer + internal probe까지 디자인하는 클래스.
	- 자동 probe 디자인 (probe_sequence=None)
	- 고정 probe_sequence 주입 모드 지원 (probe_sequence 제공)
		
	NOTE:
	- BasePrimerDesigner 쪽 opt_length, min_length, max_length, opt_gc, ... 는
	  'primer'용 조건으로 사용됩니다.
	- 이 클래스의 opt_tm, probe_min_tm, probe_max_tm, probe_* 길이/GC 조건은
	  'probe(내부 올리고)' 용으로만 사용됩니다.
	"""

	def __init__(
		self,
		template_sequence: str,
		target_start_index: int,
		target_end_index: int,
		*,
		# probe 관련
		n_probes: int = 10,
		probe_sequence: Optional[str] = None,
		probe_opt_length: int = 25,
		probe_min_length: int = 20,
		probe_max_length: int = 30,
		opt_tm: float = 60.0,		  # == probe_opt_tm (사용 예: opt_tm=probe_opt_tm_eff)
		probe_min_tm: float = 50.0,
		probe_max_tm: float = 70.0,
		probe_opt_gc: float = 45.0,
		probe_min_gc: float = 35.0,
		probe_max_gc: float = 65.0,
		probe_primer3_global_args: Optional[Dict[str, Any]] = None,
		reference_template_sequence: Optional[str] = None,
		**kwargs: Any,
	) -> None:
		# probe 개수/시퀀스
		self.n_probes: int = n_probes
		self.probe_sequence: Optional[str] = probe_sequence

		# probe 조건 (primer 조건과 별도)
		self.probe_opt_length: int = probe_opt_length
		self.probe_min_length: int = probe_min_length
		self.probe_max_length: int = probe_max_length

		self.probe_opt_tm: float = opt_tm
		self.probe_min_tm: float = probe_min_tm
		self.probe_max_tm: float = probe_max_tm

		self.probe_opt_gc: float = probe_opt_gc
		self.probe_min_gc: float = probe_min_gc
		self.probe_max_gc: float = probe_max_gc

		# 기본 primer 설정은 BasePrimerDesigner에 위임
		# - kwargs 안에 들어있는 opt_length, min_length, max_length, opt_gc, ...
		#   는 모두 primer 조건으로 사용됨
		# - 여기서 opt_tm(=probe_opt_tm)은 super()에 넘기지 않아서 primer Tm에는 영향 없음
		super().__init__(
			template_sequence=template_sequence,
			target_start_index=target_start_index,
			target_end_index=target_end_index,
			reference_template_sequence=reference_template_sequence,
			**kwargs,
		)

		# probe 사용 설정
		self.primer3_global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 1
		self.primer3_global_args["PRIMER_INTERNAL_NUM_RETURN"] = self.n_probes

		# probe 자동 디자인 / 고정 시퀀스 각각 설정
		if self.probe_sequence is not None:
			self._configure_fixed_probe()
		else:
			self._configure_probe_auto()

		# 사용자 지정 probe용 global args (primer용 primer3_global_args와 별도로 merge)
		if probe_primer3_global_args:
			self.update_primer3_global_args(probe_primer3_global_args)

	# ------------------------------------------------------------------
	# probe 자동 디자인 설정
	# ------------------------------------------------------------------
	def _configure_probe_auto(self) -> None:
		"""
		target(SNP/CpG, 1~2bp)을 반드시 포함하는 probe를 만들기 위한 설정.
		- probe start 위치는 '어떤 start를 골라도 최소 길이(probe_min_length)로 target을 항상 포함'하는 구간만 허용
		- 길이/GC/Tm 조건은 probe_* 속성을 사용
		"""
		self.update_primer3_seq_args(
				{
					'SEQUENCE_TARGET': [self.target_start_index, self.target_end_index-self.target_start_index+1],
					'SEQUENCE_INTERNAL_EXCLUDED_REGION': [[0, self.target_end_index-self.probe_min_length], [self.target_start_index+self.probe_min_length, len(self.template_sequence)-(self.target_start_index+self.probe_min_length)]]
				}
			)
		# ---- probe(내부 올리고) 특성 설정 ----
		self.update_primer3_global_args(
			{
			"PRIMER_INTERNAL_SALT_MONOVALENT": self.DEFAULT_SALT_MONOVALENT,
			"PRIMER_INTERNAL_SALT_DIVALENT": self.DEFAULT_SALT_DIVALENT,
			"PRIMER_INTERNAL_DNTP_CONC": self.DEFAULT_DNTP_CONC,
			"PRIMER_INTERNAL_DNA_CONC": self.DEFAULT_DNA_CONC,
			"PRIMER_INTERNAL_OPT_SIZE": self.probe_opt_length,
			"PRIMER_INTERNAL_MIN_SIZE": self.probe_min_length,
			"PRIMER_INTERNAL_MAX_SIZE": self.probe_max_length,
			"PRIMER_INTERNAL_OPT_TM": self.probe_opt_tm,
			"PRIMER_INTERNAL_MIN_TM": self.probe_min_tm,
			"PRIMER_INTERNAL_MAX_TM": self.probe_max_tm,
			"PRIMER_INTERNAL_OPT_GC_PERCENT": self.probe_opt_gc,
			"PRIMER_INTERNAL_MIN_GC": self.probe_min_gc,
			"PRIMER_INTERNAL_MAX_GC": self.probe_max_gc,
		}
	)


	# ------------------------------------------------------------------
	# 고정 probe 시퀀스 사용 설정
	# ------------------------------------------------------------------
	def _configure_fixed_probe(self) -> None:
		"""미리 정해진 probe_sequence를 사용하는 경우 설정."""
		probe_start = self.template_sequence.find(self.probe_sequence)
		if probe_start == -1:
			# 필요하면 여기서 ValueError로 바꿀 수도 있음
			return

		self.update_primer3_seq_args(
			{
				"SEQUENCE_INTERNAL_OLIGO": self.probe_sequence,
				# TODO: probe와 primer 사이 간격을 인자로 받도록 개선 가능
				"SEQUENCE_EXCLUDED_REGION": [
					[probe_start - 1, len(self.probe_sequence) + 2]
				],
			}
		)

		# 고정 probe의 경우 길이/Tm/GC 제약은 크게 풀어둔 상태
		# (원하면 self.probe_* 를 사용하도록 변경 가능)
		self.update_primer3_global_args(
			{
				"PRIMER_INTERNAL_SALT_MONOVALENT": self.DEFAULT_SALT_MONOVALENT,
				"PRIMER_INTERNAL_SALT_DIVALENT": self.DEFAULT_SALT_DIVALENT,
				"PRIMER_INTERNAL_DNTP_CONC": self.DEFAULT_DNTP_CONC,
				"PRIMER_INTERNAL_DNA_CONC": self.DEFAULT_DNA_CONC,
				"PRIMER_INTERNAL_MIN_SIZE": 0,
				"PRIMER_INTERNAL_MAX_SIZE": 30,
				"PRIMER_INTERNAL_MIN_TM": 0,
				"PRIMER_INTERNAL_MAX_TM": 100,
				"PRIMER_INTERNAL_MIN_GC": 0,
				"PRIMER_INTERNAL_MAX_GC": 100,
			}
		)

	def reset(self) -> None:
		super().reset()

		# probe 사용 설정 다시
		self.primer3_global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 1
		self.primer3_global_args["PRIMER_INTERNAL_NUM_RETURN"] = self.n_probes

		# probe 자동/고정 모드 다시 설정
		if self.probe_sequence is not None:
			self._configure_fixed_probe()
		else:
			self._configure_probe_auto()
	# ------------------------------------------------------------------
	# Amplicon 생성
	# ------------------------------------------------------------------
	def _build_amplicons(self) -> List[Amplicon]:
		"""
		primer + probe를 모두 Amplicon에 포함하고,
		probe가 타겟 구간 전체를 커버하는 것만 필터링.
		"""
		assert self.primer3_result is not None

		n_forward = self.primer3_result.get("PRIMER_LEFT_NUM_RETURNED", 0)
		n_reverse = self.primer3_result.get("PRIMER_RIGHT_NUM_RETURNED", 0)
		n_internal = self.primer3_result.get("PRIMER_INTERNAL_NUM_RETURNED", 0)
		n_pairs = self.primer3_result.get("PRIMER_PAIR_NUM_RETURNED", 0)

		n_designed = max(n_forward, n_reverse, n_internal, n_pairs)
		amplicons: List[Amplicon] = []

		for rank in range(n_designed):
			forward: Optional[Primer] = None
			reverse: Optional[Primer] = None
			probe: Optional[Primer] = None

			if self.primer3_result.get(f"PRIMER_LEFT_{rank}") is not None:
				forward = Primer(
					template_sequence=self.template_sequence,
					reference_template_sequence=self.reference_template_sequence,
					sequence=self.primer3_result[f"PRIMER_LEFT_{rank}_SEQUENCE"],
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					strand="forward",
					primer_type="forward",
				)

			if self.primer3_result.get(f"PRIMER_RIGHT_{rank}") is not None:
				reverse = Primer(
					template_sequence=self.template_sequence,
					reference_template_sequence=self.reference_template_sequence,
					sequence=self.primer3_result[f"PRIMER_RIGHT_{rank}_SEQUENCE"],
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					strand="reverse",
					primer_type="reverse",
				)

			if self.primer3_result.get(f"PRIMER_INTERNAL_{rank}") is not None:
				probe = Primer(
					template_sequence=self.template_sequence,
					reference_template_sequence=self.reference_template_sequence,
					sequence=self.primer3_result[f"PRIMER_INTERNAL_{rank}_SEQUENCE"],
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					strand="forward",
					primer_type="probe",
				)

			amplicon = Amplicon(
				template_sequence=self.template_sequence,
				reference_template_sequence=self.reference_template_sequence,
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				forward_primer=forward,
				reverse_primer=reverse,
				probe=probe,
			)
			# probe가 타겟 전체를 커버하는 경우만 사용 (probe 없으면 그냥 통과)
			if probe is not None:
				probe_start = probe.template_sequence.find(probe.sequence)
				probe_end = probe_start + len(probe.sequence)
				if (
					probe_start <= self.target_start_index
					and probe_end >= self.target_end_index
				):
					amplicons.append(amplicon)
			else:
				amplicons.append(amplicon)
		print(len(amplicons))
		return amplicons