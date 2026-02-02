from __future__ import annotations

from typing import Any, Dict, List, Optional

import primer3
from pcr.components import Primer, Amplicon


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
					penalty=self.primer3_result.get(f"PRIMER_LEFT_{rank}_PENALTY", 0.0),
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
					penalty=self.primer3_result.get(f"PRIMER_RIGHT_{rank}_PENALTY", 0.0),
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
# Primer + Probe
# ----------------------------------------------------------------------
class ProbePrimerDesigner(BasePrimerDesigner):
	"""
	1) probe-only 디자인 -> (forward/reverse None, probe만) Amplicon 리스트 생성
	2) 그 리스트 loop:
	   - probe를 고정(SEQUENCE_INTERNAL_OLIGO)
	   - probe 주변 gap 만큼 exclusion zone(SEQUENCE_EXCLUDED_REGION) 설정
	   - primer Tm window = probe_tm - diff 로 설정
	   - primer pair 디자인
	   - probe+primers Amplicon으로 결과 생성
	"""

	def __init__(
		self,
		template_sequence: str,
		target_start_index: int,
		target_end_index: int,
		*,
		n_probes: int = 100,
		probe_sequence: Optional[str] = None,
		probe_opt_length: int = 25,
		probe_min_length: int = 20,
		probe_max_length: int = 30,
		probe_opt_tm: float = 65.0,
		probe_min_tm: float = 60.0,
		probe_max_tm: float = 70.0,
		min_primer_probe_tm_diff: float = 5.0,
		max_primer_probe_tm_diff: float = 10.0,
		probe_opt_gc: float = 45.0,
		probe_min_gc: float = 35.0,
		probe_max_gc: float = 65.0,
		probe_primer3_global_args: Optional[Dict[str, Any]] = None,
		# ✅ NEW: probe 주변 exclusion gap (bp)
		probe_gap: int = 3,
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

		self.probe_gap = int(probe_gap)
		self._probe_primer3_global_args = probe_primer3_global_args or {}

		# ✅ summary 저장용
		self.summary: Dict[str, Any] = {
			"probe_only": {},
			"filters": {"5_g": 0, "poly_g": 0, "3_gc": 0, "target_cover_fail": 0, "template_find_fail": 0},
			"primer3_fail": {"no_pair": 0},
			"counts": {"probe_candidates": 0, "probe_after_filter": 0, "primer_runs": 0, "amplicons_final": 0},
			"qc": {},  # 네가 qc dict를 밖에서 주입하면 여기 붙이면 됨 (옵션)
		}

		super().__init__(
			template_sequence=template_sequence,
			target_start_index=target_start_index,
			target_end_index=target_end_index,
			reference_template_sequence=reference_template_sequence,
			**kwargs,
		)

	def _configure_probe_only(self) -> None:
		self.forward_primer = False
		self.reverse_primer = False
		self.primer3_global_args["PRIMER_PICK_LEFT_PRIMER"] = 0
		self.primer3_global_args["PRIMER_PICK_RIGHT_PRIMER"] = 0

		self.primer3_global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 1
		self.primer3_global_args["PRIMER_INTERNAL_NUM_RETURN"] = self.n_probes * 10
		self.primer3_global_args["PRIMER_NUM_RETURN"] = self.n_probes * 10
		print(self.n_probes * 10)
		print('xx')
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

		if self._probe_primer3_global_args:
			self.update_primer3_global_args(self._probe_primer3_global_args)

	def _build_probe_only_amplicons(self) -> List[Amplicon]:
		assert self.primer3_result is not None
		n_internal = int(self.primer3_result.get("PRIMER_INTERNAL_NUM_RETURNED", 0))

		self.summary["probe_only"] = {
			"internal_returned": n_internal,
			"internal_explain": self.primer3_result.get("PRIMER_INTERNAL_EXPLAIN"),
		}

		amplicons: List[Amplicon] = []
		for rank in range(n_internal):
			if self.primer3_result.get(f"PRIMER_INTERNAL_{rank}") is None:
				continue
			probe_seq = self.primer3_result.get(f"PRIMER_INTERNAL_{rank}_SEQUENCE")
			if not probe_seq:
				continue

			probe_penalty = self.primer3_result.get(f"PRIMER_INTERNAL_{rank}_PENALTY", 0.0)
			probe = Primer(
				template_sequence=self.template_sequence,
				reference_template_sequence=self.reference_template_sequence,
				sequence=probe_seq,
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				strand="forward",
				primer_type="probe",
				penalty=probe_penalty,
			)

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

	# ✅ NEW: probe 고정 + gap exclusion + primer tm 세팅 + primer3 실행
	def _run_primer3_with_params(
		self,
		*,
		probe: Primer,
		gap: int,
		target_tm: Optional[float],
		min_diff: float,
		max_diff: float,
	) -> Dict[str, Any]:
		# reset은 호출자가 해도 되고 여기서 해도 되는데, 여기서는 "seq/global args만" 조정한다고 가정
		# (design()에서 reset 후 호출)

		# 1) probe 고정
		self.update_primer3_seq_args({"SEQUENCE_INTERNAL_OLIGO": probe.sequence})
		self.update_primer3_global_args({"PRIMER_PICK_INTERNAL_OLIGO": 1})

		# 2) probe 주변 exclusion zone: [probe_start-gap, probe_end+gap)
		#	primer3 excluded region format: [[start, length]]
		p_start = int(getattr(probe, "binding_start_index", -1))
		if p_start < 0:
			# binding_start_index가 없다면 template에서 찾아서 계산
			p_start = self.template_sequence.find(probe.sequence)

		if p_start >= 0:
			p_len = len(probe.sequence)
			excl_start = max(0, p_start - int(gap))
			excl_end = min(len(self.template_sequence), p_start + p_len + int(gap))
			excl_len = max(0, excl_end - excl_start)
			if excl_len > 0:
				self.primer3_seq_args["SEQUENCE_EXCLUDED_REGION"] = [[excl_start, excl_len]]
			else:
				self.primer3_seq_args.pop("SEQUENCE_EXCLUDED_REGION", None)
		else:
			# 찾기 실패시 exclusion 미적용
			self.primer3_seq_args.pop("SEQUENCE_EXCLUDED_REGION", None)

		# 3) primer tm window = probe_tm - diff
		if target_tm is not None:
			self.update_primer3_global_args(
				{
					"PRIMER_OPT_TM": float(target_tm) - float(min_diff),
					"PRIMER_MIN_TM": float(target_tm) - float(max_diff),
					"PRIMER_MAX_TM": float(target_tm) - float(min_diff),
				}
			)

		# 4) run primer3
		return primer3.bindings.designPrimers(self.primer3_seq_args, self.primer3_global_args) or {}

	def _build_amplicons_from_pair_result(self, res: Dict[str, Any], probe: Primer) -> List[Amplicon]:
		n_forward = int(res.get("PRIMER_LEFT_NUM_RETURNED", 0))
		n_reverse = int(res.get("PRIMER_RIGHT_NUM_RETURNED", 0))
		n_pairs = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))

		if n_pairs == 0:
			self.summary["primer3_fail"]["no_pair"] += 1
			return []

		n_designed = max(n_forward, n_reverse, n_pairs)
		amps: List[Amplicon] = []

		for rank in range(n_designed):
			if res.get(f"PRIMER_LEFT_{rank}") is None or res.get(f"PRIMER_RIGHT_{rank}") is None:
				continue
			
			forward_penalty = res.get(f"PRIMER_LEFT_{rank}_PENALTY", 0.0)
			reverse_penalty = res.get(f"PRIMER_RIGHT_{rank}_PENALTY", 0.0)
			
			forward = Primer(
				template_sequence=self.template_sequence,
				reference_template_sequence=self.reference_template_sequence,
				sequence=res.get(f"PRIMER_LEFT_{rank}_SEQUENCE"),
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				strand="forward",
				primer_type="forward",
				penalty=forward_penalty,
			)
			reverse = Primer(
				template_sequence=self.template_sequence,
				reference_template_sequence=self.reference_template_sequence,
				sequence=res.get(f"PRIMER_RIGHT_{rank}_SEQUENCE"),
				target_start_index=self.target_start_index,
				target_end_index=self.target_end_index,
				strand="reverse",
				primer_type="reverse",
				penalty=reverse_penalty,
			)

			amps.append(
				Amplicon(
					template_sequence=self.template_sequence,
					reference_template_sequence=self.reference_template_sequence,
					target_start_index=self.target_start_index,
					target_end_index=self.target_end_index,
					forward_primer=forward,
					reverse_primer=reverse,
					probe=probe
				)
			)

		return amps

	def design(self) -> List[Amplicon]:
		final_amplicons: List[Amplicon] = []

		# 1) probe-only
		self.reset()
		self._configure_probe_only()
		self.run_primer3()
		probe_only_amplicons = self._build_probe_only_amplicons()

		self.summary["counts"]["probe_candidates"] = len(probe_only_amplicons)

		# 2) loop each probe -> primers design
		drop_counts = self.summary["filters"]
		max_poly_g = 4
		max_3prime_gc = 3

		for probe_amp in probe_only_amplicons:
			if probe_amp.probe is None:
				continue

			seq = probe_amp.probe.sequence

			# (a) target cover check
			ps = probe_amp.template_sequence.find(seq)
			if ps < 0:
				drop_counts["template_find_fail"] += 1
				continue
			pe = ps + len(seq)
			if not (ps <= self.target_start_index and pe >= self.target_end_index):
				drop_counts["target_cover_fail"] += 1
				continue

			# (b) filters (5' G, polyG, 3' GC)
			if seq.startswith("G"):
				drop_counts["5_g"] += 1
				continue
			if "G" * (max_poly_g + 1) in seq:
				drop_counts["poly_g"] += 1
				continue
			tail_5bp = seq[-5:]
			if tail_5bp.count("G") + tail_5bp.count("C") >= max_3prime_gc:
				drop_counts["3_gc"] += 1
				continue

			self.summary["counts"]["probe_after_filter"] += 1

			# probe tm
			probe_tm = getattr(probe_amp.probe, "tm", None)
			if probe_tm is not None:
				try:
					probe_tm = float(probe_tm)
				except Exception:
					probe_tm = None

			# (c) primer design mode
			self.reset()
			self.forward_primer = True
			self.reverse_primer = True
			self.primer3_global_args["PRIMER_PICK_LEFT_PRIMER"] = 1
			self.primer3_global_args["PRIMER_PICK_RIGHT_PRIMER"] = 1

			self.summary["counts"]["primer_runs"] += 1

			# ✅ run with exclusion gap + fixed probe + tm window
			res = self._run_primer3_with_params(
				probe=probe_amp.probe,
				gap=self.probe_gap,
				target_tm=probe_tm,
				min_diff=self.min_primer_probe_tm_diff,
				max_diff=self.max_primer_probe_tm_diff,
			)

			final_amplicons.extend(self._build_amplicons_from_pair_result(res, probe_amp.probe))

		self.amplicon_list = final_amplicons
		self.summary["counts"]["amplicons_final"] = len(final_amplicons)
		
		print(self.summary)

		return self.amplicon_list

	def _build_amplicons(self) -> List[Amplicon]:
		# design() override라서 기본 경로는 사용 안 함
		return []
