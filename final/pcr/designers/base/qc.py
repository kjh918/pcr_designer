from typing import List
from pcr.config.schema.app import PipelineConfig
from pcr.components.amplicon import Amplicon
from pcr.qc.thermo import ThermoChecker
from pcr.qc.blast import BlastSpecificityChecker


class BaseQCExecutor:
    def __init__(self, config: PipelineConfig):
        """
        Base QC 실행기
        - config: PipelineConfig 객체
        """
        self.config = config 
        self.checkers = []
        
        # 🔥 [Logic 추가] reference_name이 'none'이면 blast를 스킵하도록 설정
        # app.js에서 보낸 'none'은 스키마를 거쳐 config.reference_name에 담겨 있습니다.
        ref_name = getattr(self.config, "reference_name", "none")
        should_use_blast = (str(ref_name).lower() != "none")
        
        # 판별된 값을 _setup_checkers에 전달합니다.
        self._setup_checkers(blast=should_use_blast)

    def _setup_checkers(self, blast: bool = True):
        """
        체커 리스트를 초기화합니다.
        """
        # 1. ThermoChecker는 물리적 특성 검사이므로 항상 포함
        self.checkers.append(ThermoChecker(self.config.qc_criteria))
        
        # 2. blast 인자가 True이고, 유전체 정보가 있을 때만 BLAST 체커 추가
        if blast:
            print("🧬 [QC] BLAST Specificity Check is ENABLED.")
            self.checkers.append(BlastSpecificityChecker(self.config))
        else:
            print("🚀 [QC] BLAST Specificity Check is DISABLED (Skipped).")
        
    def execute(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        """
        [Core Method] 등록된 체커들을 차례로 실행하여 Amplicons를 필터링합니다.
        """
        current_candidates = amplicons
        
        for checker in self.checkers:
            if not current_candidates:
                break
            # 각 체커(ThermoChecker, BlastSpecificityChecker)의 run()을 호출
            current_candidates = checker.run(current_candidates)
            
        return current_candidates