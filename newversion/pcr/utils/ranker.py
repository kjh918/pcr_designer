from typing import List
from ..components.amplicon import Amplicon

class ProbeCentricRanker:
    def __init__(self, probe_overlap_threshold: float = 0.9):
        """
        probe_overlap_threshold: Probe 서열이 50% 이상 겹치면 같은 위치로 간주
        """
        self.threshold = probe_overlap_threshold

    def _get_total_penalty(self, amp: Amplicon) -> float:
        """
        Pair Penalty와 Internal(Probe) Penalty를 합산하여 
        세트 전체의 품질 점수를 계산
        """
        total = amp.pair_penalty  # 기본 프라이머 쌍 페널티
        
        # Probe가 있다면 Probe 페널티도 합산
        if amp.probe and hasattr(amp.probe, 'penalty'):
            total += amp.probe.penalty
            
        return total

    def _is_probe_overlapping(self, amp1: Amplicon, amp2: Amplicon) -> bool:
        """두 앰플리콘의 Probe 위치 중첩도 계산"""
        if not amp1.probe or not amp2.probe:
            return False
        
        p1, p2 = amp1.probe, amp2.probe
        intersection = max(0, min(p1.end_index, p2.end_index) - max(p1.start_index, p2.start_index))
        
        # 두 프로브 중 더 짧은 쪽의 길이를 기준으로 겹침 비율 계산
        min_len = min(p1.end_index - p1.start_index, p2.end_index - p2.start_index)
        return (intersection / min_len) >= self.threshold

    def select_diverse_probes(self, amplicons: List[Amplicon], top_k: int = 10) -> List[Amplicon]:
        """s
        1. Total Penalty(Pair + Internal) 기준으로 정렬
        2. Probe가 있는 후보를 최우선으로 선발
        3. 동일 Probe 구역 내에서는 가장 우수한 세트 1개만 남김
        """
        
        # 1. Total Penalty 계산 및 정렬 (낮을수록 우수)
        # 람다 함수 내에서 합산 점수를 기준으로 정렬합니다.
        sorted_amps = sorted(amplicons, key=self._get_total_penalty)

        # 2. 우선순위 분리 (Probe 유무)
        with_probe = [a for a in sorted_amps if a.probe is not None]

        selected: List[Amplicon] = []

        # 3. Probe 세트 선발 (위치 중복 제거)
        for amp in with_probe:
            if len(selected) >= top_k:
                break
            
            # 이미 뽑힌 세트들과 Probe 위치가 겹치는지 확인
            is_redundant = any(self._is_probe_overlapping(amp, chosen) for chosen in selected)
            
            if not is_redundant:
                # 합산 페널티 정보를 객체에 기록 (나중에 확인용)
                amp.total_penalty = self._get_total_penalty(amp)
                selected.append(amp)

        return selected