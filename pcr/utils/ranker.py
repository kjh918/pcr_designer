"""
pcr/utils/ranker.py
QC를 통과한 최종 후보군들의 우선순위를 매기고 다양성을 확보하는 모듈
"""
from typing import List
from pcr.components.amplicon import Amplicon

class ProbeCentricRanker:
    def __init__(self, probe_overlap_threshold: float = 0.9):
        """
        [MODIFIED] SNP 설계에서는 위치가 무조건 겹치므로, 
        좌표 기반의 overlap threshold는 무시하고 '서열 자체의 다름'을 기준으로 랭킹을 매깁니다.
        (파라미터는 하위 호환성을 위해 남겨둠)
        """
        pass

    def _get_total_penalty(self, amp: Amplicon) -> float:
        """Pair Penalty와 Internal(Probe) Penalty를 합산하여 세트 전체의 품질 점수를 계산"""
        total = amp.pair_penalty
        
        # Probe 페널티 안전하게 합산
        if amp.probe and getattr(amp.probe, 'penalty', 0.0):
            total += amp.probe.penalty
            
        return total

    def select_diverse_probes(self, amplicons: List[Amplicon], top_k: int = 10) -> List[Amplicon]:
        """
        1. Total Penalty 기준으로 정렬 (낮을수록 우수)
        2. Probe가 있는 후보만 선발
        3. [핵심] 완전히 동일한 서열의 Probe를 가진 세트가 이미 있다면 패스 (다양성 확보)
        """
        # 1. Total Penalty가 가장 낮은(좋은) 순서대로 정렬
        sorted_amps = sorted(amplicons, key=self._get_total_penalty)
        
        selected: List[Amplicon] = []
        seen_probe_seqs = set() # 이미 선택된 프로브의 '염기서열'을 기억하는 바구니

        for amp in sorted_amps:
            if len(selected) >= top_k:
                break
            
            if not amp.probe:
                continue
                
            probe_seq = amp.probe.sequence.upper()
            
            # 2. 프로브 서열이 바구니에 없으면(즉, 새로운 서열이면) 픽업!
            if probe_seq not in seen_probe_seqs:
                # 합산 페널티 정보를 객체에 기록 (리포팅 확인용)
                amp.total_penalty = self._get_total_penalty(amp)
                
                selected.append(amp)
                seen_probe_seqs.add(probe_seq)

        return selected