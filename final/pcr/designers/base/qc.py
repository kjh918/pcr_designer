from typing import List
from pcr.config.schema.app import PipelineConfig
from pcr.components.amplicon import Amplicon
from pcr.qc.thermo import ThermoChecker
from pcr.qc.blast import BlastSpecificityChecker


class BaseQCExecutor:
    def __init__(self, config: PipelineConfig):
        # config는 PipelineConfig 객체여야 합니다.
        self.config = config 
        self.checkers = []
        self._setup_checkers()

    def _setup_checkers(self):

        # 1. ThermoChecker는 criteria만 필요하므로 명시적으로 전달
        self.checkers.append(ThermoChecker(self.config.qc_criteria))
        
        # 2. BlastSpecificityChecker는 경로(system) 정보도 필요하므로 config 전체 전달
        self.checkers.append(BlastSpecificityChecker(self.config))
    
    def execute(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        """
        [Core Method] 자식 클래스의 super().execute()가 찾아올 목표 메서드입니다.
        """
        current_candidates = amplicons
        for checker in self.checkers:
            if not current_candidates:
                break
            # 각 체커(ThermoChecker, BlastSpecificityChecker)의 run()을 호출하여 필터링
            current_candidates = checker.run(current_candidates)
        
        return current_candidates