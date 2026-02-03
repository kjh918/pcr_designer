from typing import Optional, Any, Dict, List
from pydantic import BaseModel, Field

# =============================================================================
# 1. System Config (환경 설정 - 경로)
# =============================================================================
class QCPaths(BaseModel):
    """
    system.yaml에 정의된 외부 도구 및 DB 경로
    키 이름은 YAML 파일과 일치해야 합니다.
    """
    blast_root: Optional[str] = None
    blast_bin_dir: Optional[str] = None
    blastn: str             # 실행 파일 경로 (필수)
    blastdbcmd: Optional[str] = None
    blast_db_path: str      # BLAST DB prefix (필수)
    blast_ref_path: str     # isPCR용 Reference FASTA (필수)

    class Config:
        populate_by_name = True
        extra = "ignore" # YAML에 불필요한 키가 있어도 무시

class QCToolsConfig(BaseModel):
    """System 설정을 담는 컨테이너"""
    paths: QCPaths

# =============================================================================
# 2. Preset Config (실험 기준 - Threshold)
# =============================================================================

# (1) 공통 기준
class CommonThermoCriteria(BaseModel):
    hairpin_max_dg: float = Field(-5.0, description="Max dG for Hairpin")
    homodimer_max_dg: float = Field(-6.0, description="Max dG for Homodimer")
    heterodimer_max_dg: float = Field(-6.0, description="Max dG for Heterodimer")

class CommonSpecificityCriteria(BaseModel):
    blast_max_alignments: int = Field(50, alias="max_alignments")
    blast_identity_threshold: float = Field(90.0, alias="min_identity")
    blast_word_size: int = Field(12, alias="word_size")
    min_amp_size: int = Field(50, alias="min_amp_len")
    max_amp_size: int = Field(3000, alias="max_amp_len")

    class Config:
        populate_by_name = True

# (2) Primer/Probe/AS-PCR 개별 기준
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
    
    class Config:
        populate_by_name = True

# (3) 통합 Criteria 모델 (Root)
class QCCriteria(BaseModel):
    """default.yaml 매핑"""
    # Thermo
    hairpin_min_dg: float = -5.0
    homodimer_min_dg: float = -6.0
    heterodimer_min_dg: float = -6.0
    
    # Specificity
    min_identity: float = 95.0
    min_hit_length: int = 15
    blast_max_alignments: int = 1000
    min_amp_size: int = 50
    max_amp_size: int = 3000
    
    # Sub-criteria
    primer: PrimerQCCriteria = Field(default_factory=PrimerQCCriteria)
    probe: ProbeQCCriteria = Field(default_factory=ProbeQCCriteria)
    as_pcr: Optional[ASPCRQCCriteria] = None
    
    # IsPCR 사용 여부
    use_ispcr_check: bool = True

    class Config:
        extra = "ignore"

# =============================================================================
# 3. QC 결과 및 상태 모델
# =============================================================================
class QCResult(BaseModel):
    """개별 QC 항목의 통과 여부 및 상세 정보"""
    passed: bool
    value: Any = None
    threshold: Any = None
    reason: Optional[str] = None
    data: Optional[Dict[str, Any]] = None

class AmpliconQCStatus(BaseModel):
    """
    모든 QC 결과를 통합하여 저장하는 컨테이너.
    QCExecutor가 이 구조체에 결과를 담아 Amplicon.qc_status에 저장합니다.
    """
    thermo: Optional[Dict[str, Any]] = None      # ThermoChecker 결과
    specificity: Optional[Dict[str, Any]] = None # BlastSpecificityChecker 결과
    ispcr: Optional[Dict[str, Any]] = None       # IsPcrChecker 결과
    overall_passed: bool = False                 # 최종 통과 여부

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
    product_size: int
    fwd_hit: Optional[BlastHit] = None
    rev_hit: Optional[BlastHit] = None
    is_target: bool = False