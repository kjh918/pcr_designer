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

    # pcr/qc/executor.py (BaseQCExecutor)

    def _setup_checkers(self):
        # 1. ThermoChecker는 항상 실행 (물리적 특성 체크)
        self.checkers.append(ThermoChecker(self.config.qc_criteria))
        
        # 2. BLAST DB 경로가 설정되어 있을 때만 BlastSpecificityChecker 추가
        # PipelineConfig의 qc_criteria 혹은 system 설정에 blast_db가 있는지 확인
        blast_db = getattr(self.config.qc_criteria, "blast_db", None)
        
        if blast_db:
            self.checkers.append(BlastSpecificityChecker(self.config))
        else:
            # 로그를 남겨서 사용자가 BLAST가 생략되었음을 알게 합니다.
            import logging
            logging.getLogger(__name__).info("🚀 BLAST DB path is None. Skipping BlastSpecificityChecker.")
        
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