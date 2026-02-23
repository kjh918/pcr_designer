"""
pcr/factory.py
assay_type을 받아 적절한 Designer와 QC 파이프라인을 생성하고 실행하는 통합 진입점(Facade).
Design -> QC -> Ranking 순서로 파이프라인을 관장합니다.
"""
from typing import Optional, List, Dict, Any

# Config & Schemas
from .config.loader import load_pipeline_config
from .config.schema.app import PipelineConfig
from .designers.base.schema import BaseDesignInput, BaseDesignOutput

# QPCR
from .designers.qpcr.designer import QPCRPrimerDesigner
from .designers.qpcr.schema import QPCRDesignInput
from .designers.qpcr.qc import QPCRQCExecutor

# # AS-PCR
# from .designers.as_pcr.designer import ASPCRPrimerDesigner
# from .designers.as_pcr.schema import ASPCRDesignInput
# from .designers.as_pcr.qc import ASPCRQCExecutor

# # MS-PCR
# from .designers.ms_pcr.designer import MSPCRPrimerDesigner
# from .designers.ms_pcr.schema import MSPCRDesignInput
# from .designers.ms_pcr.qc import MSPCRQCExecutor

# Ranker
from .utils.ranker import ProbeCentricRanker

# ─────────────────────────────────────────────────────────
# 통합 Registry: Assay 타입에 따른 (Input, Designer, QC) 세트 매핑
# ─────────────────────────────────────────────────────────
PIPELINE_MAP = {
    "qpcr":   (QPCRDesignInput, QPCRPrimerDesigner, QPCRQCExecutor),
#     "as_pcr": (ASPCRDesignInput, ASPCRPrimerDesigner, ASPCRQCExecutor),
#     "ms_pcr": (MSPCRDesignInput, MSPCRPrimerDesigner, MSPCRQCExecutor),
}
class PCRFactory:
    def __init__(self, config: Optional[PipelineConfig] = None):
        self.config = config

    def run(
        self,
        assay_type: str,
        name: str,
        target_start: int,
        target_end: int,
        reference_name: str = "hg38",
        template_sequence: Optional[str] = None,
        reference_sequence: Optional[str] = None,
        top_k: int = 5,
        run_qc: bool = True,
        overrides: dict = None,
        **extra_input_kwargs,
    ) -> BaseDesignOutput:
        
        assay_type = assay_type.lower()
        if assay_type not in PIPELINE_MAP:
            raise ValueError(f"Unknown assay_type: '{assay_type}'. Choose from {list(PIPELINE_MAP.keys())}")

        InputClass, DesignClass, QCClass = PIPELINE_MAP[assay_type]

        if not self.config:
            self.config = load_pipeline_config(
                base_yaml_path="pcr/config/base_pcr.yaml",
                system_yaml_path="pcr/config/system.yaml",
                assay_type=assay_type,
                user_overrides=overrides
            )

        # ---------------------------------------------------------
        # 1. Input 생성
        # ---------------------------------------------------------
        design_input = InputClass(
            name=name,
            template_sequence=template_sequence,
            target_start=target_start,
            target_end=target_end,
            config=self.config,
            reference_name=reference_name,
            reference_sequence=reference_sequence,
            overrides=overrides or {},
            **extra_input_kwargs,
        )

        # ---------------------------------------------------------
        # 2. Design 실행
        # ---------------------------------------------------------
        designer = DesignClass(design_input)
        output = designer.design()

        if output.status != "success" or not output.amplicons:
            return output
            
        initial_count = len(output.amplicons)

        # ---------------------------------------------------------
        # 3. QC 실행 및 탈락 사유 추적
        # ---------------------------------------------------------
        if run_qc:
            qc_executor = QCClass(self.config)
            
            # [수정] QC 실행
            qc_passed_amplicons = qc_executor.execute(output.amplicons)
            
            # [핵심 추가] QC 통계 수집
            # 가정: qc_executor 내부에 self.qc_stats = {"blast_fail": 0, "thermo_fail": 0, ...} 와 같은 딕셔너리가 존재함
            qc_stats = getattr(qc_executor, "qc_stats", {})
            
            # 터미널에 즉시 출력 (디버깅 용도)
            print(f"\n🔍 [QC STATS] Total Analyzed: {initial_count} | Passed: {len(qc_passed_amplicons)}")
            if qc_stats:
                for reason, count in qc_stats.items():
                    print(f"   => Failed due to '{reason}': {count}")
            
            # JSON 결과 로그에 추가 (웹 출력용)
            output.log_messages.append(f"QC Passed: {len(qc_passed_amplicons)} / {initial_count}")
            if qc_stats:
                output.log_messages.append(f"QC Fail Reasons: {qc_stats}")
                
        else:
            qc_passed_amplicons = output.amplicons
            output.log_messages.append("QC Skipped.")

        # ---------------------------------------------------------
        # 4. Ranking
        # ---------------------------------------------------------
        if qc_passed_amplicons:
            ranker = ProbeCentricRanker(probe_overlap_threshold=0.9)
            final_amplicons = ranker.select_diverse_probes(qc_passed_amplicons, top_k=top_k)
            output.amplicons = final_amplicons
            output.log_messages.append(f"Ranker Selected: {len(final_amplicons)} (Top K: {top_k})")
        else:
            output.amplicons = []
            output.status = "fail"
            output.error_msg = "All amplicons failed QC."
            # 모두 탈락했을 때 상세 사유 로그 추가
            if run_qc and getattr(qc_executor, "qc_stats", {}):
                output.error_msg += f" Details: {qc_executor.qc_stats}"

        return output