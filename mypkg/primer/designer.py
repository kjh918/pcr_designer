from __future__ import annotations

from typing import Any, Dict, List, Optional

import primer3
from primer.pcr_components import Primer, Amplicon

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
    - 자동 probe 디자인 (probe=True)
    - 혹은 고정 probe_sequence 주입 모드 지원
    """

    def __init__(
        self,
        template_sequence: str,
        target_start_index: int,
        target_end_index: int,
        *,
        n_probes: int = 10,
        probe_sequence: Optional[str] = None,
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
        **kwargs: Any,
    ) -> None:
        self.n_probes: int = n_probes
        self.probe_sequence: Optional[str] = probe_sequence

        # 기본 primer 설정은 그대로 상속
        super().__init__(
            template_sequence=template_sequence,
            target_start_index=target_start_index,
            target_end_index=target_end_index,
            opt_length=opt_length,
            min_length=min_length,
            max_length=max_length,
            opt_tm=opt_tm,
            min_tm=min_tm,
            max_tm=max_tm,
            opt_gc=opt_gc,
            min_gc=min_gc,
            max_gc=max_gc,
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

    def _configure_probe_auto(self) -> None:
        """
        타겟 전체를 커버하는 probe 자동 디자인을 위한 primer3 설정.
        """
        target_len = self.target_end_index - self.target_start_index + 1

        self.update_primer3_seq_args(
            {
                "SEQUENCE_TARGET": [self.target_start_index, target_len],
                "SEQUENCE_INTERNAL_EXCLUDED_REGION": [
                    [0, self.target_end_index - self.max_length],
                    [
                        self.target_start_index + self.max_length,
                        len(self.template_sequence)
                        - (self.target_start_index + self.max_length),
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
                "PRIMER_INTERNAL_OPT_SIZE": self.opt_length,
                "PRIMER_INTERNAL_MIN_SIZE": self.min_length,
                "PRIMER_INTERNAL_MAX_SIZE": self.max_length,
                "PRIMER_INTERNAL_OPT_TM": self.opt_tm,
                "PRIMER_INTERNAL_MIN_TM": self.min_tm,
                "PRIMER_INTERNAL_MAX_TM": self.max_tm,
                "PRIMER_INTERNAL_OPT_GC_PERCENT": self.opt_gc,
                "PRIMER_INTERNAL_MIN_GC": self.min_gc,
                "PRIMER_INTERNAL_MAX_GC": self.max_gc,
            }
        )

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
                "SEQUENCE_EXCLUDED_REGION": [[probe_start - 1, len(self.probe_sequence) + 2]],
            }
        )

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

        return amplicons
