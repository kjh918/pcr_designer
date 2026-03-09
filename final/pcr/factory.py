"""
pcr/factory.py
assay_type을 받아 적절한 Designer와 QC 파이프라인을 생성하고 실행하는 통합 진입점(Facade).
Design -> QC -> Ranking 순서로 파이프라인을 관장합니다.

[Update]
- PipelineConfig의 get_reference() 메서드를 사용하여 Genome 경로(FASTA, BLAST DB)를 안전하게 조회합니다.
- 조회된 경로를 QC 파라미터 등에 주입(Injection)합니다.
"""
from typing import Optional, List, Dict, Any
import logging

# Config & Schemas
from .config.loader import load_pipeline_config
from .config.schema.app import PipelineConfig
from .designers.base.schema import BaseDesignInput, BaseDesignOutput

# QPCR
from .designers.qpcr.designer import QPCRPrimerDesigner
from .designers.qpcr.schema import QPCRDesignInput
from .designers.qpcr.qc import QPCRQCExecutor

 # AS-PCR (추후 활성화)
from .designers.as_pcr.designer import ASPCRPrimerDesigner
from .designers.as_pcr.schema import ASPCRDesignInput
from .designers.as_pcr.qc import ASPCRQCExecutor

# # MS-PCR (추후 활성화)
from .designers.ms_pcr.designer import MSPCRPrimerDesigner
from .designers.ms_pcr.schema import MSPCRDesignInput
from .designers.ms_pcr.qc import MSPCRQCExecutor

# Ranker
from .utils.ranker import ProbeCentricRanker

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────
# 통합 Registry: Assay 타입에 따른 (Input, Designer, QC) 세트 매핑
# ─────────────────────────────────────────────────────────
PIPELINE_MAP = {
    "qpcr":   (QPCRDesignInput, QPCRPrimerDesigner, QPCRQCExecutor),
    "aspcr": (ASPCRDesignInput, ASPCRPrimerDesigner, ASPCRQCExecutor),
    "mspcr": (MSPCRDesignInput, MSPCRPrimerDesigner, MSPCRQCExecutor),
}

class PCRFactory:
    def __init__(self, config: Optional[PipelineConfig] = None):
        self.config = config

    def _resolve_genome_paths(self, reference_name: str):
        """
        [FIXED] Config의 get_reference 메서드를 사용하여 경로를 안전하게 가져옵니다.
        ReferenceConfig 객체에서 FASTA 및 BLAST DB 경로를 추출하여 파이프라인 설정에 주입합니다.
        """
        if not self.config:
            logger.warning("Config is not loaded yet.")
            return

        try:
            # ✅ 1. PipelineConfig에 정의된 헬퍼 메서드 사용 (AttributeError 해결)
            ref_config = self.config.get_reference(reference_name)
        except ValueError as e:
            # system.yaml에 해당 reference가 정의되지 않은 경우
            logger.warning(f"Genome path resolution failed: {e}")
            return

        # ✅ 2. Pydantic 모델 속성 접근
        fasta_path = ref_config.fasta_path
        blast_db_path = getattr(ref_config, "blast_db", None)

        logger.info(f"Resolved paths for {reference_name}: FASTA={fasta_path}, BLAST={blast_db_path}")

        # 3. QC Criteria에 BLAST DB 경로 주입 (QCExecutor가 사용)
        if blast_db_path:
            # config.qc_criteria가 Pydantic 모델인 경우
            if hasattr(self.config.qc_criteria, "blast_db"):
                self.config.qc_criteria.blast_db = blast_db_path
            # 혹시 dict인 경우 대비
            elif isinstance(self.config.qc_criteria, dict):
                self.config.qc_criteria["blast_db"] = blast_db_path
            
            # (옵션) QC Criteria에 genome 이름도 명시
            if hasattr(self.config.qc_criteria, "genome_assembly"):
                 self.config.qc_criteria.genome_assembly = reference_name

        # 4. System 설정에 현재 사용할 FASTA 경로 업데이트 (Designer가 참조할 경우)
        # self.config.system이 SystemConfig 객체라고 가정
        if fasta_path and hasattr(self.config, "system"):
             if hasattr(self.config.system, "current_fasta"):
                 self.config.system.current_fasta = fasta_path

    def run(
        self,
        assay_type: str,
        name: str,
        target_start: int,
        target_end: int,
        reference_name: str = "hg38",
        template_sequence: Optional[str] = None,
        reference_sequence: Optional[str] = None,
        template_genomic_start: Optional[int] = 0,
        template_genomic_end: Optional[int] = 0, 
        top_k: int = 5,
        run_qc: bool = True,
        overrides: dict = None,
        **extra_input_kwargs,
    ) -> BaseDesignOutput:
        
        assay_type = assay_type.lower()
        if assay_type not in PIPELINE_MAP:
            raise ValueError(f"Unknown assay_type: '{assay_type}'. Choose from {list(PIPELINE_MAP.keys())}")

        InputClass, DesignClass, QCClass = PIPELINE_MAP[assay_type]

        # 1. Config 로드 (없을 경우)
        if not self.config:
            # 주의: load_pipeline_config가 반환하는 객체 구조가 PipelineConfig와 일치해야 함
            self.config = load_pipeline_config(
                base_yaml_path="pcr/config/base_pcr.yaml",
                system_yaml_path="pcr/config/system.yaml", 
                assay_type=assay_type,
                user_overrides=overrides
            )

        # 2. [핵심] Genome 경로 해석 및 Config 업데이트
        # 사용자 입력(hg38) -> 실제 경로(/storage/...)로 변환하여 Config에 심어줌
        self._resolve_genome_paths(reference_name)

        # 3. Input 객체 생성
        design_input = InputClass(
            name=name,
            template_sequence=template_sequence,
            target_start=target_start,
            target_end=target_end,
            config=self.config,
            reference_name=reference_name,
            reference_sequence=reference_sequence,
            template_genomic_start=template_genomic_start, # 💡 객체에 직접 꽂아줍니다!
            template_genomic_end=template_genomic_end, # 💡 객체에 직접 꽂아줍니다!
            overrides=overrides or {},
            **extra_input_kwargs,
        )
        # 4. Design 실행
        designer = DesignClass(design_input)
        output = designer.design()
        if output.status != "success" or not output.amplicons:
            return output
            
        initial_count = len(output.amplicons)

        # 5. QC 실행
        if run_qc:
            # Config에 이미 _resolve_genome_paths를 통해 blast_db가 주입된 상태임
            qc_executor = QCClass(self.config)
            qc_passed_amplicons = qc_executor.execute(output.amplicons)
            
            # QC 통계 수집 (Executor 내부에 stats가 있다고 가정)
            qc_stats = getattr(qc_executor, "qc_stats", {})
            
            # 터미널 출력 (디버깅용)
            print(f"\n🔍 [QC STATS] Genome: {reference_name} | Total: {initial_count} -> Passed: {len(qc_passed_amplicons)}")
            if qc_stats:
                for reason, count in qc_stats.items():
                    print(f"   => Failed due to '{reason}': {count}")
            
            output.log_messages.append(f"QC Passed: {len(qc_passed_amplicons)} / {initial_count}")
            if qc_stats:
                output.log_messages.append(f"QC Fail Reasons: {qc_stats}")
                
        else:
            qc_passed_amplicons = output.amplicons
            output.log_messages.append("QC Skipped.")
            
        # 6. Ranking (🔥 AS-PCR 모드 우회 처리 추가)
        if qc_passed_amplicons:
            if assay_type == "aspcr":
                output.amplicons = qc_passed_amplicons
                output.log_messages.append(f"AS-PCR Mode: Bypassed Probe Ranker. Total {len(qc_passed_amplicons)} amplicons kept.")
            else:
                ranker = ProbeCentricRanker(probe_overlap_threshold=0.9)
                final_amplicons = ranker.select_diverse_probes(qc_passed_amplicons, top_k=top_k)
                output.amplicons = final_amplicons
                output.log_messages.append(f"Ranker Selected: {len(final_amplicons)} (Top K: {top_k})")
        else:
            output.amplicons = []
            output.status = "fail"
            output.error_msg = "All amplicons failed QC."
            if run_qc and getattr(qc_executor, "qc_stats", {}):
                output.error_msg += f" Details: {qc_executor.qc_stats}"

        return output