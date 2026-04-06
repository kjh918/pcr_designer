"""
pcr/factory.py
assay_type을 받아 적절한 Designer와 QC 파이프라인을 생성하고 실행하는 통합 진입점(Facade).
Design -> QC -> Ranking 순서로 파이프라인을 관장합니다.

[Update]
- PipelineConfig의 get_reference() 메서드를 사용하여 Genome 경로(FASTA, BLAST DB)를 안전하게 조회합니다.
- 조회된 경로를 QC 파라미터 등에 주입(Injection)합니다.
- 🔥 Factory에서 공통으로 파싱 규칙([REF, ALT])을 적용하여 분리한 후, 각 Designer의 Input Schema로 데이터를 주입합니다.
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

# AS-PCR
from .designers.as_pcr.designer import ASPCRPrimerDesigner
from .designers.as_pcr.schema import ASPCRDesignInput
from .designers.as_pcr.qc import ASPCRQCExecutor

# MS-PCR
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
        Config의 get_reference 메서드를 사용하여 경로를 안전하게 가져옵니다.
        ReferenceConfig 객체에서 FASTA 및 BLAST DB 경로를 추출하여 파이프라인 설정에 주입합니다.
        """
        if not self.config:
            logger.warning("Config is not loaded yet.")
            return

        try:
            # PipelineConfig에 정의된 헬퍼 메서드 사용
            ref_config = self.config.get_reference(reference_name)
        except ValueError as e:
            logger.warning(f"Genome path resolution failed: {e}")
            return

        fasta_path = ref_config.fasta_path
        blast_db_path = getattr(ref_config, "blast_db", None)
        logger.info(f"Resolved paths for {reference_name}: FASTA={fasta_path}, BLAST={blast_db_path}")

        # QC Criteria에 BLAST DB 경로 주입
        if blast_db_path:
            if hasattr(self.config.qc_criteria, "blast_db"):
                self.config.qc_criteria.blast_db = blast_db_path
            elif isinstance(self.config.qc_criteria, dict):
                self.config.qc_criteria["blast_db"] = blast_db_path
            
            if hasattr(self.config.qc_criteria, "genome_assembly"):
                 self.config.qc_criteria.genome_assembly = reference_name

        # System 설정에 현재 사용할 FASTA 경로 업데이트
        if fasta_path and hasattr(self.config, "system"):
             if hasattr(self.config.system, "current_fasta"):
                 self.config.system.current_fasta = fasta_path

    def run(
        self,
        assay_type: str,
        name: str,
        target_start: Optional[int] = None, # 파싱을 통해 자동 획득하므로 Optional
        target_end: Optional[int] = None,   # 파싱을 통해 자동 획득하므로 Optional
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
            self.config = load_pipeline_config(
                base_yaml_path="pcr/config/base_pcr.yaml",
                system_yaml_path="pcr/config/system.yaml", 
                assay_type=assay_type,
                user_overrides=overrides
            )

        # 2. Genome 경로 해석 및 Config 업데이트
        self._resolve_genome_paths(reference_name)

        # 🔥 3. 모든 Designer 공통: 대괄호 파싱 및 REF/ALT 분리
        # 괄호가 존재하면 매핑된 DesignClass의 정적 메서드를 호출하여 해석합니다.
        if template_sequence and '[' in template_sequence and ']' in template_sequence:
            template_seq, reference_seq, target_indices = DesignClass.parse_sequence_with_brackets(template_sequence)
            template_sequence = template_seq
            # 파싱된 레퍼런스 서열 자동 주입 (사용자가 별도로 지정하지 않은 경우)
            if not reference_sequence:
                reference_sequence = reference_seq
            
            # 파싱된 타겟의 양 끝 위치를 할당 (0-based)
            if target_indices:
                target_start = min(target_indices) - 1
                target_end = max(target_indices)
            
            # 각 Designer의 Input Schema에서 알아서 가공할 수 있도록 원본 인덱스도 통째로 넘겨줌
            extra_input_kwargs["target_indices"] = target_indices
            extra_input_kwargs["target_cpg_indices"] = target_indices # MS-PCR 하위 호환성 유지

        # 방어 로직: 그래도 타겟 좌표가 없다면 에러
        if target_start is None or target_end is None:
            raise ValueError("Target positions (target_start, target_end) are required, or sequence must contain brackets [ ].")
        
        
        # 4. Input 객체 생성 (각 설계 기법별 Schema로 데이터가 흘러 들어갑니다)
        design_input = InputClass(
            name=name,
            template_sequence=template_sequence,
            target_start=target_start,
            target_end=target_end,
            config=self.config,
            reference_name=reference_name,
            reference_sequence=reference_sequence,
            template_genomic_start=template_genomic_start,
            template_genomic_end=template_genomic_end,
            overrides=overrides or {},
            **extra_input_kwargs,
        )

        # 5. Design 실행
        print(design_input)
        designer = DesignClass(design_input)
        output = designer.design()

        if output.status != "success" or not output.amplicons:
            return output
            
        initial_count = len(output.amplicons)

        # 6. QC 실행
        if run_qc:
            qc_executor = QCClass(self.config, reference_name=reference_name)
            qc_passed_amplicons = qc_executor.execute(output.amplicons)
            
            qc_stats = getattr(qc_executor, "qc_stats", {})
            
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
            
        # 7. Ranking
        if qc_passed_amplicons:
            # MS-PCR, AS-PCR은 세트나 페널티 구조가 달라서 랭커를 우회
            if assay_type in ["aspcr", "mspcr",'qpcr']:
                output.amplicons = qc_passed_amplicons
                output.log_messages.append(f"AS/MS-PCR Mode: Bypassed Probe Ranker. Total {len(qc_passed_amplicons)} amplicons kept.")
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