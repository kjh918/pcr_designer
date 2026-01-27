from pydantic import BaseModel, Field, ConfigDict
from typing import Optional

# --------------------------------------------------------------------------
# 1. 하위 컴포넌트 정의 (YAML의 섹션들과 1:1 매핑)
# --------------------------------------------------------------------------

class QCPaths(BaseModel):
    """
    YAML: qc_params -> paths 섹션
    """
    BLAST_ROOT: str
    BLAST_BIN_DIR: str
    BLASTN: str
    BLASTDBCMD: str
    BLASTN_DB: str
    BLASTN_REF: str

class BaseQCCriteria(BaseModel):
    """공통 QC 기준"""
    # YAML 키(대문자) -> Python 변수(소문자) alias
    hairpin_min_dg: float = Field(-5.0, alias="HAIRPIN_MIN_DG")
    homodimer_min_dg: float = Field(-6.0, alias="HOMODIMER_MIN_DG")
    heterodimer_min_dg: float = Field(-6.0, alias="HETERODIMER_MIN_DG")
    
    blast_identity_threshold: float = Field(90.0, alias="BLAST_IDENTITY_THRESHOLD")
    blast_length_threshold: int = Field(13, alias="BLAST_LENGTH_THRESHOLD")
    blast_max_alignments: int = Field(50, alias="BLAST_MAX_ALIGNMENTS")
    
    min_amp_bp: int = Field(50, alias="MIN_AMP_BP")
    max_amp_bp: int = Field(300, alias="MAX_AMP_BP")

    model_config = ConfigDict(populate_by_name=True)

class PrimerQCCriteria(BaseQCCriteria):
    """YAML: qc_params -> primer_criteria 섹션"""
    primer_max_diff_tm: float = Field(3.0, alias="PRIMER_MAX_DIFF_TM")
    primer_min_diff_tm: float = Field(0.0, alias="PRIMER_MIN_DIFF_TM")

class ProbeQCCriteria(BaseQCCriteria):
    """YAML: qc_params -> probe_criteria 섹션"""
    probe_max_diff_tm: float = Field(8.0, alias="PROBE_MAX_DIFF_TM")
    probe_min_diff_tm: float = Field(6.0, alias="PROBE_MIN_DIFF_TM")
    avoid_5_prime_g: bool = Field(True, alias="AVOID_5_PRIME_G")

class ASPCRQCCriteria(BaseQCCriteria):
    """YAML: qc_params -> as_pcr_criteria 섹션"""
    check_3_prime_specificity: bool = Field(True, alias="CHECK_3_PRIME_SPECIFICITY")

# --------------------------------------------------------------------------
# 2. 통합 QC Config (구조 변경됨)
# --------------------------------------------------------------------------

class QCParams(BaseModel):
    """
    [FIXED] YAML의 계층 구조를 반영하도록 수정됨.
    Flat fields 대신 Nested Model을 사용합니다.
    """
    # YAML에 'paths'라는 키가 있으면 자동으로 QCPaths 모델로 매핑됨
    paths: QCPaths  
    
    # YAML에 '*_criteria' 키가 있으면 매핑됨 (Optional 처리로 유연성 확보)
    primer_criteria: Optional[PrimerQCCriteria] = None
    probe_criteria: Optional[ProbeQCCriteria] = None
    as_pcr_criteria: Optional[ASPCRQCCriteria] = None

    # Helper 메서드: 이제 self.paths로 바로 접근 가능하므로 단순 리턴
    def get_paths(self) -> QCPaths:
        return self.paths

    def get_primer_criteria(self) -> PrimerQCCriteria:
        if not self.primer_criteria:
            return PrimerQCCriteria() # 기본값 반환
        return self.primer_criteria

    def get_probe_criteria(self) -> ProbeQCCriteria:
        if not self.probe_criteria:
            return ProbeQCCriteria()
        return self.probe_criteria