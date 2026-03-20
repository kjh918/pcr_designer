from typing import Optional, Any, Dict, List
from pydantic import BaseModel, Field

# =============================================================================
# 1. System Config (환경 설정 - 경로)
# =============================================================================
class QCPaths(BaseModel):
    blast_root: Optional[str] = None
    blast_bin_dir: Optional[str] = None
    blastn: str          # 실행 파일 경로 (필수)
    blastdbcmd: Optional[str] = None
    blast_db_path: str    # BLAST DB prefix (필수)
    blast_ref_path: str  # isPCR용 Reference FASTA (필수)

    model_config = {
        "populate_by_name": True,
        "extra": "ignore"
    }

class QCToolsConfig(BaseModel):
    paths: QCPaths

# =============================================================================
# 2. Preset Config (실험 기준 - Threshold)
# =============================================================================
class CommonThermoCriteria(BaseModel):
    hairpin_max_dg: float = Field(-5.0, description="Max dG for Hairpin")
    homodimer_max_dg: float = Field(-6.0, description="Max dG for Homodimer")
    heterodimer_max_dg: float = Field(-6.0, description="Max dG for Heterodimer")

class CommonSpecificityCriteria(BaseModel):
    blast_max_alignments: int = Field(50, alias="max_alignments")
    blast_identity_threshold: float = Field(90.0, alias="min_identity")
    blast_word_size: int = Field(12, alias="word_size")
    min_amp_size: int = Field(50, alias="min_amp_len")
    max_amp_size: int = Field(1000, alias="max_amp_len")

    model_config = {"populate_by_name": True}

class PrimerQCCriteria(BaseModel):
    max_diff_tm: float = Field(3.0)
    min_diff_tm: float = 0.0

class ProbeQCCriteria(BaseModel):
    min_primer_probe_tm_diff: float = Field(5.0)
    max_primer_probe_tm_diff: float = 10.0
    max_probe_poly_g: int = Field(3)
    max_probe_3_end_gc: int = Field(2)
    avoid_5_prime_g: bool = Field(True)

class ASPCRQCCriteria(BaseModel):
    check_3_prime_specificity: bool = Field(True, alias="CHECK_3_PRIME_SPECIFICITY")
        
    model_config = {"populate_by_name": True}

class QCCriteria(BaseModel):
    """default.yaml 매핑 (Root)"""
    hairpin_min_dg: float = Field(-5.0)
    homodimer_min_dg: float = Field(-6.0)
    heterodimer_min_dg: float = Field(-6.0)
        
    min_identity: float = Field(95.0, alias="blast_identity_threshold")
    min_hit_length: int = Field(15)
    blast_max_alignments: int = Field(1000, alias="max_alignments")
    min_amp_size: int = Field(50, alias="min_amp_len")
    max_amp_size: int = Field(3000, alias="max_amp_len")
    
    min_query_coverage: float = Field(0.8, alias="min_query_coverage")
    min_identity_threshold: float = Field(80.0, alias="min_id_threshold")
    end_match_tolerance: int = Field(1, alias="end_match_tolerance")
        
    primer: PrimerQCCriteria = Field(default_factory=PrimerQCCriteria)
    probe: ProbeQCCriteria = Field(default_factory=ProbeQCCriteria)
    as_pcr: Optional[ASPCRQCCriteria] = None
        
    model_config = {"populate_by_name": True, "extra": "ignore"}
        
# =============================================================================
# 3. 🔥 QC 결과 및 상태 모델 (Workflow 확장 지원)
# =============================================================================
class QCModuleResult(BaseModel):
    """
    개별 QC 모듈(체커)의 실행 결과를 규격화한 모델입니다.
    예: module_name="thermo", is_pass=False, messages=["High Hairpin"]
    """
    module_name: str
    is_pass: bool
    messages: List[str] = Field(default_factory=list)
    metrics: Dict[str, Any] = Field(default_factory=dict) # 상세 수치 (dG 값 등)

class AmpliconQCStatus(BaseModel):
    """
    Amplicon 객체 내부에 부착되어 QC 워크플로우 전체의 상태를 추적합니다.
    어떤 새로운 QC 모듈이 추가되더라도 이 구조 하나로 모두 수용 가능합니다.
    """
    is_pass: bool = True
    fail_reasons: List[str] = Field(default_factory=list)
    
    # 각 체커별 상세 결과를 딕셔너리로 보관 (key: module_name)
    modules: Dict[str, QCModuleResult] = Field(default_factory=dict)

    def add_result(self, module_name: str, is_pass: bool, messages: List[str] = None, metrics: Dict[str, Any] = None):
        """
        [상태 업데이트 메서드]
        체커(Checker)가 이 메서드를 호출하여 자신의 결과를 등록하면,
        자동으로 전체 is_pass 상태와 fail_reasons가 갱신됩니다.
        """
        msgs = messages or []
        
        # 1. 모듈별 결과 생성 및 저장
        result = QCModuleResult(
            module_name=module_name,
            is_pass=is_pass,
            messages=msgs,
            metrics=metrics or {}
        )
        self.modules[module_name] = result
        
        # 2. 글로벌 상태(Workflow) 자동 갱신
        if not is_pass:
            self.is_pass = False
            self.fail_reasons.extend(msgs)


# =============================================================================
# 4. BLAST 관련 모델
# =============================================================================
class BlastHit(BaseModel):
    qseqid: str
    sseqid: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    qseq: str
    sseq: str

    @property
    def strand(self) -> str:
        return "+" if self.sstart < self.send else "-"
    @property
    def genomic_start(self) -> int:
        return min(self.sstart, self.send)
    @property
    def genomic_end(self) -> int:
        return max(self.sstart, self.send)

class OffTargetAmplicon(BaseModel):
    chrom: str
    start: int
    end: int
    reference_sequence: Optional[str] = None
    template_sequence: Optional[str] = None
    product_size: int
    fwd_hit: Optional[BlastHit] = None
    rev_hit: Optional[BlastHit] = None
    is_target: bool = False
    probe_binds: bool = False