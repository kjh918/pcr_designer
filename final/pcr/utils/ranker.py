"""
pcr/utils/ranker.py
QC를 통과한 최종 후보군들의 우선순위를 매기고 다양성을 확보하는 모듈
"""
from typing import List
from pcr.components.amplicon import Amplicon

class ProbeCentricRanker:
    def __init__(self, probe_overlap_threshold: float = 0.9):
        """
        probe_overlap_threshold: Probe 서열이 90% 이상 겹치면 같은 위치로 간주
        """
        self.threshold = probe_overlap_threshold

    def _get_total_penalty(self, amp: Amplicon) -> float:
        """Pair Penalty와 Internal(Probe) Penalty를 합산하여 세트 전체의 품질 점수를 계산"""
        total = amp.pair_penalty
        
        # [MODIFIED] Probe 페널티 안전하게 합산
        if amp.probe and getattr(amp.probe, 'penalty', 0.0):
            total += amp.probe.penalty
            
        return total

    def _is_probe_overlapping(self, amp1: Amplicon, amp2: Amplicon) -> bool:
        """두 앰플리콘의 Probe 위치 중첩도 계산"""
        if not amp1.probe or not amp2.probe:
            return False
        
        p1, p2 = amp1.probe, amp2.probe
        
        # [MODIFIED] end_index가 속성으로 없을 경우를 대비한 동적 계산
        p1_end = getattr(p1, 'end_index', p1.start_index + len(p1.sequence))
        p2_end = getattr(p2, 'end_index', p2.start_index + len(p2.sequence))
        
        intersection = max(0, min(p1_end, p2_end) - max(p1.start_index, p2.start_index))
        
        min_len = min(p1_end - p1.start_index, p2_end - p2.start_index)
        if min_len == 0: 
            return False
            
        return (intersection / min_len) >= self.threshold

    def select_diverse_probes(self, amplicons: List[Amplicon], top_k: int = 10) -> List[Amplicon]:
        """
        1. Total Penalty 기준으로 정렬
        2. Probe가 있는 후보 최우선 선발
        3. 동일 구역(Overlap) 내에서는 가장 우수한 세트 1개만 남김
        """
        sorted_amps = sorted(amplicons, key=self._get_total_penalty)
        with_probe = [a for a in sorted_amps if a.probe is not None]

        selected: List[Amplicon] = []

        for amp in with_probe:
            if len(selected) >= top_k:
                break
            
            is_redundant = any(self._is_probe_overlapping(amp, chosen) for chosen in selected)
            
            if not is_redundant:
                # 합산 페널티 정보를 객체에 기록 (리포팅 확인용)
                amp.total_penalty = self._get_total_penalty(amp)
                selected.append(amp)

        return selected