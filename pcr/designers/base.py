from __future__ import annotations

from typing import Any, Dict, List, Optional

import primer3
from pcr.components import Primer, Amplicon
from pcr.config import abc

# ----------------------------------------------------------------------
# Base
# ----------------------------------------------------------------------
class BasePrimerDesigner:
	DEFAULT_SALT_MONOVALENT: float = 50.0
	DEFAULT_SALT_DIVALENT: float = 1.5
	DEFAULT_DNTP_CONC: float = 0.6
	DEFAULT_DNA_CONC: float = 50.0

	def _init_primer3_args(self) -> None:
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
		self.template_sequence: str = template_sequence
		self.reference_template_sequence: str = reference_template_sequence or template_sequence
		self.target_start_index: int = int(target_start_index)
		self.target_end_index: int = int(target_end_index)

		self.min_amplicon_length: int = int(min_amplicon_length)
		self.max_amplicon_length: int = int(max_amplicon_length)
		self.n_primers: int = int(n_primers)
		self.max_tm_difference: float = float(max_tm_difference)

		self.forward_primer: bool = bool(forward_primer)
		self.reverse_primer: bool = bool(reverse_primer)

		self.opt_length: int = int(opt_length)
		self.min_length: int = int(min_length)
		self.max_length: int = int(max_length)

		self.opt_tm: float = float(opt_tm)
		self.min_tm: float = float(min_tm)
		self.max_tm: float = float(max_tm)

		self.opt_gc: float = float(opt_gc)
		self.min_gc: float = float(min_gc)
		self.max_gc: float = float(max_gc)

		self.primer3_seq_args: Dict[str, Any] = {
			"SEQUENCE_ID": "PRIMER",
			"SEQUENCE_TEMPLATE": self.template_sequence,
		}
		self.primer3_global_args: Dict[str, Any] = {
			"PRIMER_TASK": "generic",
			"PRIMER_NUM_RETURN": self.n_primers,
			"PRIMER_PICK_LEFT_PRIMER": int(self.forward_primer),
			"PRIMER_PICK_RIGHT_PRIMER": int(self.reverse_primer),
			"PRIMER_PICK_INTERNAL_OLIGO": 0,
		}

		self._configure_primer_common()

		if primer3_seq_args:
			self.update_primer3_seq_args(primer3_seq_args)
		if primer3_global_args:
			self.update_primer3_global_args(primer3_global_args)

		self.primer3_result: Optional[Dict[str, Any]] = None
		self.amplicon_list: List[Amplicon] = []

	def _configure_primer_common(self) -> None:
		target_len = self.target_end_index - self.target_start_index + 1
		self.update_primer3_seq_args({"SEQUENCE_TARGET": [self.target_start_index, target_len]})

		# ✅ primer3 expects [[min,max]]
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
				"PRIMER_PRODUCT_SIZE_RANGE": [[self.min_amplicon_length, self.max_amplicon_length]],
			}
		)

	def update_primer3_seq_args(self, args: Dict[str, Any]) -> None:
		self.primer3_seq_args.update(args)

	def update_primer3_global_args(self, args: Dict[str, Any]) -> None:
		self.primer3_global_args.update(args)

	def run_primer3(self) -> None:
		self.primer3_result = primer3.bindings.designPrimers(
			seq_args=self.primer3_seq_args,
			global_args=self.primer3_global_args,
		)
		self.amplicon_list = self._build_amplicons()

	def _build_amplicons(self) -> List[Amplicon]:
		raise NotImplementedError

	def design(self) -> List[Amplicon]:
		self.run_primer3()
		return self.amplicon_list

	def reset(self) -> None:
		self._init_primer3_args()
		self.primer3_result = None
		self.amplicon_list = []


# ----------------------------------------------------------------------
# Primer only
# ----------------------------------------------------------------------
class PrimerDesigner(BasePrimerDesigner):
	def _build_amplicons(self) -> List[Amplicon]:
		assert self.primer3_result is not None

		n_forward = int(self.primer3_result.get("PRIMER_LEFT_NUM_RETURNED", 0))
		n_reverse = int(self.primer3_result.get("PRIMER_RIGHT_NUM_RETURNED", 0))
		n_pairs = int(self.primer3_result.get("PRIMER_PAIR_NUM_RETURNED", 0))

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

			amplicons.append(
				Amplicon(
					template_sequence=self.template_sequence,
					reference_template_sequence=self.reference_template_sequence,
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					forward_primer=forward,
					reverse_primer=reverse,
					probe=None,
				)
			)

		return amplicons


# ----------------------------------------------------------------------
# Primer + Probe (✅ 네가 말한 방식 그대로)
# ----------------------------------------------------------------------
class ProbePrimerDesigner(BasePrimerDesigner):
	"""
	Flow (네가 말한 그대로):

	1) probe-only 디자인 (forward/reverse = False)
	   -> forward/reverse 가 None인 Amplicon(probe만) 리스트를 만든다
	2) 그 Amplicon 리스트를 loop:
	   - 각 probe를 고정(SEQUENCE_INTERNAL_OLIGO)
	   - 그 probe에 맞춰 primer 조건을 세팅(예: primer Tm = probe_tm - diff)
	   - primer pair 디자인
	   - 결과 Amplicon에 (probe + forward + reverse)로 저장
	"""

	def __init__(
		self,
		template_sequence: str,
		target_start_index: int,
		target_end_index: int,
		*,
		n_probes: int = 10,
		probe_sequence: Optional[str] = None,
		probe_opt_length: int = 25,
		probe_min_length: int = 20,
		probe_max_length: int = 30,
		probe_opt_tm: float = 60.0,
		probe_min_tm: float = 50.0,
		probe_max_tm: float = 70.0,
		min_primer_probe_tm_diff: float = 6.0,
		max_primer_probe_tm_diff: float = 8.0,
		probe_opt_gc: float = 45.0,
		probe_min_gc: float = 35.0,
		probe_max_gc: float = 65.0,
		probe_primer3_global_args: Optional[Dict[str, Any]] = None,
		reference_template_sequence: Optional[str] = None,
		**kwargs: Any,
	) -> None:
		self.n_probes = int(n_probes)
		self.probe_sequence = probe_sequence

		self.probe_opt_length = int(probe_opt_length)
		self.probe_min_length = int(probe_min_length)
		self.probe_max_length = int(probe_max_length)

		self.probe_opt_tm = float(probe_opt_tm)
		self.probe_min_tm = float(probe_min_tm)
		self.probe_max_tm = float(probe_max_tm)

		self.min_primer_probe_tm_diff = float(min_primer_probe_tm_diff)
		self.max_primer_probe_tm_diff = float(max_primer_probe_tm_diff)

		self.probe_opt_gc = float(probe_opt_gc)
		self.probe_min_gc = float(probe_min_gc)
		self.probe_max_gc = float(probe_max_gc)

		self._probe_primer3_global_args = probe_primer3_global_args or {}

		super().__init__(
			template_sequence=template_sequence,
			target_start_index=target_start_index,
			target_end_index=target_end_index,
			reference_template_sequence=reference_template_sequence,
			**kwargs,
		)

	# ------------------------------------------------------------------
	# probe-only primer3 설정
	# ------------------------------------------------------------------
	def _configure_probe_only(self) -> None:
		# ✅ primer는 끄고 probe만 켠다
		self.forward_primer = False
		self.reverse_primer = False
		self.primer3_global_args["PRIMER_PICK_LEFT_PRIMER"] = 0
		self.primer3_global_args["PRIMER_PICK_RIGHT_PRIMER"] = 0

		self.primer3_global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 1
		self.primer3_global_args["PRIMER_INTERNAL_NUM_RETURN"] = self.n_probes
		self.primer3_global_args["PRIMER_NUM_RETURN"] = self.n_probes  # explain/limit용

		# target 포함 조건 + excluded region (기존 네 로직 유지)
		self.update_primer3_seq_args(
			{
				"SEQUENCE_TARGET": [
					self.target_start_index,
					self.target_end_index - self.target_start_index + 1,
				],
				"SEQUENCE_INTERNAL_EXCLUDED_REGION": [
					[0, self.target_end_index - self.probe_min_length],
					[
						self.target_start_index + self.probe_min_length,
						len(self.template_sequence) - (self.target_start_index + self.probe_min_length),
					],
				],
			}
		)

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

		# 사용자 override
		if self._probe_primer3_global_args:
			self.update_primer3_global_args(self._probe_primer3_global_args)

	# ------------------------------------------------------------------
	# 고정 probe 시퀀스(선택)
	# ------------------------------------------------------------------
	def _configure_fixed_probe(self) -> None:
		if not self.probe_sequence:
			return
		self.update_primer3_seq_args({"SEQUENCE_INTERNAL_OLIGO": self.probe_sequence})
		self.primer3_global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 1

	# ------------------------------------------------------------------
	# probe-only 결과를 "probe만 있는 Amplicon"으로 만든다 (forward/reverse None)
	# ------------------------------------------------------------------
	def _build_probe_only_amplicons(self) -> List[Amplicon]:
		assert self.primer3_result is not None
		n_internal = int(self.primer3_result.get("PRIMER_INTERNAL_NUM_RETURNED", 0))
		amplicons: List[Amplicon] = []

		for rank in range(n_internal):
			if self.primer3_result.get(f"PRIMER_INTERNAL_{rank}") is None:
				continue
			probe_seq = self.primer3_result.get(f"PRIMER_INTERNAL_{rank}_SEQUENCE")
			if not probe_seq:
				continue

			probe = Primer(
				template_sequence=self.template_sequence,
				reference_template_sequence=self.reference_template_sequence,
				sequence=probe_seq,
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				strand="forward",
				primer_type="probe",
			)

			# ✅ forward/reverse None 인 amplicon 생성
			amplicons.append(
				Amplicon(
					template_sequence=self.template_sequence,
					reference_template_sequence=self.reference_template_sequence,
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					forward_primer=None,
					reverse_primer=None,
					probe=probe,
				)
			)

		return amplicons

	# ------------------------------------------------------------------
	# probe 고정 후 primer pair 결과를 probe+primer Amplicon으로 만든다
	# ------------------------------------------------------------------
	def _build_amplicons_from_pair_result(self, res: Dict[str, Any], probe: Primer) -> List[Amplicon]:
		n_forward = int(res.get("PRIMER_LEFT_NUM_RETURNED", 0))
		n_reverse = int(res.get("PRIMER_RIGHT_NUM_RETURNED", 0))
		n_pairs = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))

		n_designed = max(n_forward, n_reverse, n_pairs)
		amps: List[Amplicon] = []

		for rank in range(n_designed):
			if res.get(f"PRIMER_LEFT_{rank}") is None or res.get(f"PRIMER_RIGHT_{rank}") is None:
				continue

			forward = Primer(
				template_sequence=self.template_sequence,
				reference_template_sequence=self.reference_template_sequence,
				sequence=res.get(f"PRIMER_LEFT_{rank}_SEQUENCE"),
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				strand="forward",
				primer_type="forward",
			)
			reverse = Primer(
				template_sequence=self.template_sequence,
				reference_template_sequence=self.reference_template_sequence,
				sequence=res.get(f"PRIMER_RIGHT_{rank}_SEQUENCE"),
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				strand="reverse",
				primer_type="reverse",
			)

			amps.append(
				Amplicon(
					template_sequence=self.template_sequence,
					reference_template_sequence=self.reference_template_sequence,
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					forward_primer=forward,
					reverse_primer=reverse,
					probe=probe,
				)
			)

		return amps

	# ------------------------------------------------------------------
	# ✅ 핵심: 네가 말한 전체 플로우를 design()에서 수행
	# ------------------------------------------------------------------
	def design(self) -> List[Amplicon]:
		final_amplicons: List[Amplicon] = []

		# (A) probe_sequence가 있으면: 그 probe로만 primer 디자인
		if self.probe_sequence:
			# primer 디자인 모드로 초기화
			self.reset()
			self.forward_primer = True
			self.reverse_primer = True
			self.primer3_global_args["PRIMER_PICK_LEFT_PRIMER"] = 1
			self.primer3_global_args["PRIMER_PICK_RIGHT_PRIMER"] = 1
			self._configure_fixed_probe()

			# probe Primer 객체 만들고
			probe_obj = Primer(
				template_sequence=self.template_sequence,
				reference_template_sequence=self.reference_template_sequence,
				sequence=self.probe_sequence,
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				strand="forward",
				primer_type="probe",
			)

			# primer3 실행
			res = primer3.bindings.designPrimers(self.primer3_seq_args, self.primer3_global_args) or {}
			final_amplicons = self._build_amplicons_from_pair_result(res, probe_obj)

			self.amplicon_list = final_amplicons
			return self.amplicon_list

		# (B) 1) probe-only 디자인 실행 -> probe만 있는 amplicon 리스트 생성
		self.reset()
		self._configure_probe_only()
		self.run_primer3()
		probe_only_amplicons = self._build_probe_only_amplicons()
		
		# (C) 2) probe-only amplicon loop -> 각 probe 고정 후 primer pair 디자인
		for probe_amp in probe_only_amplicons:
			if probe_amp.probe is None:
				continue

			# probe가 target 커버하는 것만 (안전)
			ps = probe_amp.template_sequence.find(probe_amp.probe.sequence)
			if ps < 0:
				continue
			pe = ps + len(probe_amp.probe.sequence)
			if not (ps <= self.target_start_index and pe >= self.target_end_index):
				continue

			# probe tm 기반 primer Tm window 조절 (네가 쓰던 방식)
			probe_tm = getattr(probe_amp.probe, "tm", None)
			if probe_tm is not None:
				try:
					probe_tm = float(probe_tm)
				except Exception:
					probe_tm = None

			# primer 디자인 모드로 reset
			self.reset()
			self.forward_primer = True
			self.reverse_primer = True
			self.primer3_global_args["PRIMER_PICK_LEFT_PRIMER"] = 1
			self.primer3_global_args["PRIMER_PICK_RIGHT_PRIMER"] = 1

			# probe 고정 주입
			self.update_primer3_seq_args({"SEQUENCE_INTERNAL_OLIGO": probe_amp.probe.sequence})
			self.update_primer3_global_args({"PRIMER_PICK_INTERNAL_OLIGO": 1})

			# ✅ primer Tm = probe_tm - diff (네 qPCRdesigner 로직)
			if probe_tm is not None:
				self.update_primer3_global_args(
					{
						"PRIMER_OPT_TM": probe_tm - self.min_primer_probe_tm_diff,
						"PRIMER_MIN_TM": probe_tm - self.max_primer_probe_tm_diff,
						"PRIMER_MAX_TM": probe_tm - self.min_primer_probe_tm_diff,
					}
				)

			res = primer3.bindings.designPrimers(self.primer3_seq_args, self.primer3_global_args) or {}
			final_amplicons.extend(self._build_amplicons_from_pair_result(res, probe_amp.probe))

		self.amplicon_list = final_amplicons
		return self.amplicon_list

	# Base가 호출하는 hook은 여기서는 사용 안 하므로 안전장치
	def _build_amplicons(self) -> List[Amplicon]:
		assert self.primer3_result is not None
		# design()을 override했기 때문에 기본경로로는 안 들어오는 게 정상
		return []
