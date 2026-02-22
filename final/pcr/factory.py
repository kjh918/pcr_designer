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
        """
        config가 주어지지 않으면 기본 yaml 경로에서 로드합니다.
        (실무에서는 외부에서 assay_type에 맞게 load_pipeline_config()로 생성한 객체를 주입받는 것을 권장)
        """
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
        """
        [MODIFIED] 파이프라인 오케스트레이션: Design -> QC -> Rank
        """
        assay_type = assay_type.lower()
        if assay_type not in PIPELINE_MAP:
            raise ValueError(f"Unknown assay_type: '{assay_type}'. Choose from {list(PIPELINE_MAP.keys())}")

        InputClass, DesignClass, QCClass = PIPELINE_MAP[assay_type]

        # Config Lazy Loading (주입되지 않았을 경우)
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
        # 3. QC 실행 (Batch 최적화 적용)
        # ---------------------------------------------------------
        if run_qc:
            # Assay 특화 QC 파이프라인 초기화
            qc_executor = QCClass(self.config)
            
            # [MODIFIED] for문으로 하나씩 돌리지 않고, 리스트 전체를 넘겨 배치 처리(BLAST 등) 수행
            qc_passed_amplicons = qc_executor.execute(output.amplicons)
            
            output.log_messages.append(f"QC Passed: {len(qc_passed_amplicons)} / {initial_count}")
        else:
            qc_passed_amplicons = output.amplicons
            output.log_messages.append("QC Skipped.")

        # ---------------------------------------------------------
        # 4. Ranking (QC를 통과한 결과물 대상)
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

        return output