# pcr/designers/as_pcr.py
from __future__ import annotations

from typing import Any, Dict, List, Optional

from pcr.designers.base import BasePrimerDesigner
from pcr.components import Primer, Amplicon


def _build_left_anchor_excluded_region(
    *,
    template_len: int,
    target_index: int,
    min_len: int,
    max_len: int,
) -> List[List[int]]:
    """
    LEFT primer의 3' end == target_index가 되도록
    LEFT primer start 가능한 범위만 남기고 나머지 start를 제외하는 excluded region 생성.

    LEFT primer는 forward이므로:
      start = end - L + 1
      end를 target_index로 고정하면,
      start ∈ [target_index - max_len + 1, target_index - min_len + 1]
    """
    t = int(target_index)
    Lmin = int(min_len)
    Lmax = int(max_len)

    if template_len <= 0:
        return []

    start_min = max(0, t - Lmax + 1)
    start_max = max(0, t - Lmin + 1)

    if start_min > start_max:
        return [[0, template_len]]

    excluded: List[List[int]] = []

    # 앞쪽 제외: [0, start_min)
    if start_min > 0:
        excluded.append([0, start_min])

    # 뒤쪽 제외: (start_max, end]
    tail_start = start_max + 1
    if tail_start < template_len:
        excluded.append([tail_start, template_len - tail_start])

    return excluded


def _build_right_anchor_excluded_region(
    *,
    template_len: int,
    target_index: int,
    min_len: int,
    max_len: int,
) -> List[List[int]]:
    """
    RIGHT primer의 3' end == target_index가 되도록 start 범위를 제한.

    primer3에서 PRIMER_RIGHT의 location(start,len)에서 start는 "leftmost index"로 주어짐.
    reverse primer의 3' end가 target_index가 되려면, (template 상에서)
      binding_end == target_index
      binding_end = start + len - 1
      => start = target_index - (len - 1)

    len이 [min_len, max_len] 이므로 start 범위는
      start ∈ [target_index - (max_len - 1), target_index - (min_len - 1)]
           = [target_index - max_len + 1, target_index - min_len + 1]

    즉 LEFT 앵커와 start 범위는 동일하게 나온다.
    """
    # 결과적으로 LEFT와 동일 범위가 된다(길이 가변이면).
    # 다만 의미는 "RIGHT primer의 end 고정"이다.
    return _build_left_anchor_excluded_region(
        template_len=template_len,
        target_index=target_index,
        min_len=min_len,
        max_len=max_len,
    )


class AsPcrDesigner(BasePrimerDesigner):
    """
    (1) LEFT anchored run: forward primer 3' end == target_index
    (2) RIGHT anchored run: reverse primer 3' end == target_index

    입력:
      - target_start_index/target_end_index: BasePrimerDesigner용 타겟(1bp면 start=end)
      - target_index: 앵커링 기준 위치(보통 start와 동일하게 넘김)
      - ref_allele/alt_allele: 이후 WT/ALT 생성 단계에서 사용(지금 단계에서는 저장만)
    """

    def __init__(
        self,
        *,
        template_sequence: str,
        reference_template_sequence: str,
        target_start_index: int,
        target_end_index: int,
        target_index: int,
        ref_allele: str,
        alt_allele: str,
        min_amplicon_length: int,
        max_amplicon_length: int,
        n_primers: int = 50,
        primer3_global_args: Optional[Dict[str, Any]] = None,
    ) -> None:
        super().__init__(
            template_sequence=template_sequence,
            reference_template_sequence=reference_template_sequence,
            target_start_index=int(target_start_index),
            target_end_index=int(target_end_index),
            min_amplicon_length=int(min_amplicon_length),
            max_amplicon_length=int(max_amplicon_length),
            n_primers=int(n_primers),
            forward_primer=True,
            reverse_primer=True,
            primer3_global_args=primer3_global_args,
        )

        self.target_index = int(target_index)
        self.ref_allele = (ref_allele or "").strip().upper()
        self.alt_allele = (alt_allele or "").strip().upper()

        if len(self.ref_allele) != 1 or len(self.alt_allele) != 1:
            raise ValueError("ref_allele/alt_allele must be 1bp each (e.g. A/G).")

    # ---------------------------
    # Anchor 모드 설정
    # ---------------------------
    def _apply_left_anchor(self) -> None:
        """
        LEFT primer의 3' end가 target_index가 되도록 탐색을 제한.
        ✅ 여기서는 'LEFT 앵커'만 강제하고 RIGHT는 자유롭게 둔다.
        """
        excluded = _build_left_anchor_excluded_region(
            template_len=len(self.template_sequence),
            target_index=self.target_index,
            min_len=self.min_length,
            max_len=self.max_length,
        )

        # primer3는 SEQUENCE_EXCLUDED_REGION이 전체 primer 탐색에 영향을 준다.
        # 따라서 "LEFT만" 완벽히 제한하는 옵션은 없고,
        # 실전적으로는 2-pass + 사후검증으로 안정화한다.
        self.update_primer3_seq_args({"SEQUENCE_EXCLUDED_REGION": excluded})

        # pair는 계속 뽑되, left쪽이 앵커 범위 밖이면 나오기 어렵게 된다.
        self.update_primer3_global_args(
            {
                "PRIMER_PICK_LEFT_PRIMER": 1,
                "PRIMER_PICK_RIGHT_PRIMER": 1,
                "PRIMER_PICK_INTERNAL_OLIGO": 0,
            }
        )

    def _apply_right_anchor(self) -> None:
        """
        RIGHT primer의 3' end가 target_index가 되도록 탐색 제한.
        """
        excluded = _build_right_anchor_excluded_region(
            template_len=len(self.template_sequence),
            target_index=self.target_index,
            min_len=self.min_length,
            max_len=self.max_length,
        )
        self.update_primer3_seq_args({"SEQUENCE_EXCLUDED_REGION": excluded})

        self.update_primer3_global_args(
            {
                "PRIMER_PICK_LEFT_PRIMER": 1,
                "PRIMER_PICK_RIGHT_PRIMER": 1,
                "PRIMER_PICK_INTERNAL_OLIGO": 0,
            }
        )

    # ---------------------------
    # primer3 결과 → Amplicon
    # ---------------------------
    def _build_amplicons(self) -> List[Amplicon]:
        assert self.primer3_result is not None
        res = self.primer3_result

        amplicons: List[Amplicon] = []
        n_pairs = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0) or 0)

        for i in range(n_pairs):
            f_loc = res.get(f"PRIMER_LEFT_{i}")
            r_loc = res.get(f"PRIMER_RIGHT_{i}")
            f_seq = res.get(f"PRIMER_LEFT_{i}_SEQUENCE")
            r_seq = res.get(f"PRIMER_RIGHT_{i}_SEQUENCE")

            if not f_loc or not r_loc or not f_seq or not r_seq:
                continue

            f_start, f_len = int(f_loc[0]), int(f_loc[1])
            r_start, r_len = int(r_loc[0]), int(r_loc[1])

            forward = Primer(
                template_sequence=self.template_sequence,
                reference_template_sequence=self.reference_template_sequence,
                sequence=f_seq,
                strand="forward",
                primer_type="forward",
                target_start_index=self.target_start_index,
                target_end_index=self.target_end_index,
                binding_start_index=f_start,
                binding_end_index=f_start + f_len - 1,
            )

            reverse = Primer(
                template_sequence=self.template_sequence,
                reference_template_sequence=self.reference_template_sequence,
                sequence=r_seq,
                strand="reverse",
                primer_type="reverse",
                target_start_index=self.target_start_index,
                target_end_index=self.target_end_index,
                binding_start_index=r_start,
                binding_end_index=r_start + r_len - 1,
            )

            amplicons.append(
                Amplicon(
                    template_sequence=self.template_sequence,
                    reference_template_sequence=self.reference_template_sequence,
                    target_start_index=self.target_start_index,
                    target_end_index=self.target_end_index,
                    forward_primer=forward,
                    reverse_primer=reverse,
                )
            )

        return amplicons

    # ---------------------------
    # 2-pass design
    # ---------------------------
    def design(self) -> Dict[str, List[Amplicon]]:
        results: Dict[str, List[Amplicon]] = {}

        # LEFT anchored
        self.reset()
        self._apply_left_anchor()
        self.run_primer3()
        results["LEFT_ANCHORED"] = self._build_amplicons()

        # RIGHT anchored
        self.reset()
        self._apply_right_anchor()
        self.run_primer3()
        results["RIGHT_ANCHORED"] = self._build_amplicons()

        return results
