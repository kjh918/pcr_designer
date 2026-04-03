"""
pcr/designers/base/qc.py
모든 PCR 기법(qPCR, AS-PCR, MS-PCR)에서 공통으로 사용하는 QC 실행기 기저 클래스입니다.
유전체 참조 이름이 'none'인 경우 BLAST 특이성 검사를 자동으로 생략합니다.
"""
from typing import List
from pcr.config.schema.app import PipelineConfig
from pcr.components.amplicon import Amplicon
from pcr.qc.thermo import ThermoChecker
from pcr.qc.blast import BlastSpecificityChecker

class BaseQCExecutor:
    """
    프라이머의 물리적 특성(열역학) 및 서열 특이성(BLAST)을 검증하는 기본 클래스입니다.
    """
    def __init__(self, config: PipelineConfig, reference_name: str = "none"):
        """
        QC 실행기를 초기화하고 설정에 따라 필요한 체커들을 등록합니다.
        """
        self.config = config 
        self.checkers = []
        
        should_use_blast = (str(reference_name).lower() != "none")
        self._setup_checkers(blast=should_use_blast)

    def _setup_checkers(self, blast: bool = True):
        """
        사용될 QC 모듈을 등록합니다. 
        자식 클래스에서 이 메서드를 오버라이딩할 때 반드시 'blast' 인자를 포함해야 합니다.
        """
        # 1. ThermoChecker (Hairpin, Dimer 등)는 항상 실행합니다.
        self.checkers.append(ThermoChecker(self.config.qc_criteria))
        
        # 2. BLAST 체커는 유전체 정보가 있고 사용자가 선택한 경우에만 추가합니다.
        if blast:
            print("🧬 [QC] BLAST Specificity Check is ENABLED.")
            self.checkers.append(BlastSpecificityChecker(self.config))
        else:
            print("🚀 [QC] BLAST Specificity Check is DISABLED (Skipped).")
        
    def execute(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        """
        [Core Method] 등록된 모든 체커를 순차적으로 실행하여 앰플리콘 리스트를 필터링합니다.
        """
        current_candidates = amplicons
        
        for checker in self.checkers:
            if not current_candidates:
                break
            # 각 체커의 run()을 호출하여 결과(is_qc_pass 등)를 업데이트합니다.
            current_candidates = checker.run(current_candidates)
            
        return current_candidates