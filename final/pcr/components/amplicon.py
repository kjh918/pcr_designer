from __future__ import annotations
from typing import Any, Dict, List, Optional, Union
from pydantic import BaseModel, Field, model_validator, ConfigDict
import primer3

# 내부 모듈 임포트
from .region import GenomicRegion, SequenceChange
from .primer import Primer, Probe

class Amplicon(BaseModel):
    """
    PCR Amplicon 객체 (Pydantic Model)
    - Forward/Reverse Primer로 정의되는 증폭 산물
    - 변이(Variation) 분석 및 물리적 성질(Tm, Size) 자동 계산
    """
    # -------------------------------------------------------------------------
    # 1. 필수 입력 필드
    # -------------------------------------------------------------------------
    id: str = Field(..., description="Amplicon ID")
    forward: Primer
    reverse: Primer
    set_id: Optional[str] = None
    probe: Optional[Probe] = None
        
    template_sequence: str = Field(..., description="전체 템플릿 서열 (Context 포함)")
    target_start_index: int = Field(-1, description="타겟 영역 시작 (0-based)")
    target_end_index: int = Field(-1, description="타겟 영역 끝 (0-based)")
        
    reference_sequence: str = Field("", description="변이 분석을 위한 Reference 서열")
    reference_id: str = Field("", description="Chromosome or Gene ID")
        
    # -------------------------------------------------------------------------
    # 2. 결과 및 상태 필드 (QC 결과 저장 필드 추가)
    # -------------------------------------------------------------------------
    pair_penalty: float = 0.0
    total_penalty: float = 0.0
    allele_type: Optional[str] = None
    is_qc_pass: bool = False
    qc_log: str = ""
    off_target_count: int = 0
        
    # ✅ [ADDED/FIXED] QC 결과를 담기 위한 핵심 필드들
    # ThermoChecker가 dG 값 등을 저장하는 곳
    thermo_stats: Dict[str, Any] = Field(default_factory=dict)
    # BlastSpecificityChecker가 In-silico PCR 결과를 저장하는 곳
    blast_stats: Dict[str, Any] = Field(default_factory=dict)
    
    # 하위 호환성을 위한 범용 QC 상세 데이터
    qc_details: Dict[str, Any] = Field(default_factory=dict)
        
    # QCExecutor가 AmpliconQCStatus 객체를 저장하는 곳
    qc_status: Optional[Any] = None 
        
    # 내부 계산 필드
    sequence: str = ""        # 실제 증폭된 Amplicon 서열
    ref_sequence_clip: str = "" # Amplicon 위치에 해당하는 Reference 서열
    product_size: int = 0
    tm: float = 0.0
    gc: float = 0.0
    region: Optional[GenomicRegion] = None
    genomic_pos: str = ""       # Amplicon 게놈 좌표 (chr:start-end)
    alignment_visual: List[str] = Field(default_factory=list, description="웹 UI 렌더링용 Alignment 다이어그램")    
    # 변이 분석 결과
    all_changes: List[SequenceChange] = Field(default_factory=list)
    target_change_count: int = 0
    mismatch_count: int = 0

    # Pydantic 설정: 임의의 타입 허용
    model_config = ConfigDict(arbitrary_types_allowed=True)

    # -------------------------------------------------------------------------
    # 3. 초기화 로직
    # -------------------------------------------------------------------------
    @model_validator(mode='after')
    def compute_properties(self) -> 'Amplicon':
        """객체 생성 직후 물리적 성질 및 변이 분석 수행"""
        self._calculate_properties()
        self.reanalyze_variations()
        return self

    def _calculate_properties(self):
        """기본 물성(Tm, Size) 및 좌표 계산"""
        if self.forward.start_index is not None and self.reverse.end_index is not None:
            # 1. 길이 계산
            self.product_size = self.reverse.end_index - self.forward.start_index
            
            # 2. 서열 추출
            if self.template_sequence:
                self.sequence = self.template_sequence[self.forward.start_index : self.reverse.end_index]
                try:
                    self.tm = primer3.calc_tm(self.sequence, mv_conc=50, dv_conc=1.5, dntp_conc=0.6, dna_conc=50)
                    g_count = self.sequence.upper().count('G')
                    c_count = self.sequence.upper().count('C')
                    self.gc = round(((g_count + c_count) / len(self.sequence)) * 100, 2)
                except Exception:
                    self.tm = 0.0
                    self.gc = 0.0

            # 3. Reference 서열 추출
            if self.reference_sequence and len(self.reference_sequence) >= len(self.template_sequence):
                self.ref_sequence_clip = self.reference_sequence[self.forward.start_index : self.reverse.end_index]

        # 4. Genomic Region 매핑
        if self.forward.region and self.reverse.region:
            f_reg = self.forward.region
            r_reg = self.reverse.region
            coords = [f_reg.start, f_reg.end, r_reg.start, r_reg.end]
            self.region = GenomicRegion(
                chrom=f_reg.chrom, 
                start=min(coords), 
                end=max(coords), 
                strand=f_reg.strand
            )

    def _is_cpg_context(self, index: int, ref_base: str) -> bool:
        """CpG Site 판별"""
        if not self.ref_sequence_clip: return False
        seq_len = len(self.ref_sequence_clip)
        
        if ref_base == 'C':
            if index + 1 < seq_len:
                return self.ref_sequence_clip[index + 1].upper() == 'G'
        elif ref_base == 'G':
            if index > 0:
                return self.ref_sequence_clip[index - 1].upper() == 'C'
        return False

    def reanalyze_variations(self):
        """변이(Mismatch) 재분석 로직"""
        self.all_changes = []
        self.target_change_count = 0
        self.mismatch_count = 0

        if not self.sequence or not self.ref_sequence_clip:
            return
        
        if len(self.sequence) != len(self.ref_sequence_clip):
            return

        fwd_len = len(self.forward.sequence) if self.forward.sequence else 0
        rev_len = len(self.reverse.sequence) if self.reverse.sequence else 0
        total_len = len(self.sequence)

        # 타겟 상대 좌표 계산
        amp_start_abs = self.forward.start_index
        rel_t_start = -1
        rel_t_end = -1

        if self.target_start_index != -1 and amp_start_abs is not None:
            rel_t_start = self.target_start_index - amp_start_abs
            rel_t_end = self.target_end_index - amp_start_abs

        for i, (alt, ref) in enumerate(zip(self.sequence, self.ref_sequence_clip)):
            if alt.upper() != ref.upper():
                ref_u, alt_u = ref.upper(), alt.upper()
                
                on_target = False
                if rel_t_start != -1:
                    if rel_t_start <= i < rel_t_end:
                        on_target = True

                r_type = "internal"
                if i < fwd_len:
                    r_type = "forward_primer"
                elif i >= (total_len - rev_len):
                    r_type = "reverse_primer"

                is_conv = False
                if (ref_u == 'C' and alt_u == 'T') or (ref_u == 'G' and alt_u == 'A'):
                    if self._is_cpg_context(i, ref_u):
                        is_conv = True

                change = SequenceChange(
                    position=i,
                    ref_base=ref_u,
                    alt_base=alt_u,
                    region_type=r_type,
                    on_target=on_target,
                    is_bisulfite_conversion=is_conv
                )
                self.all_changes.append(change)

                if on_target:
                    self.target_change_count += 1
                else:
                    self.mismatch_count += 1

    # -------------------------------------------------------------------------
    # API Methods
    # -------------------------------------------------------------------------
    def set_target_region(self, start: int, end: int):
        self.target_start_index = start
        self.target_end_index = end
        self.reanalyze_variations()

    def manually_approve_change(self, position: int):
        for change in self.all_changes:
            if change.position == position:
                if not change.on_target:
                    change.on_target = True
                    self.mismatch_count = max(0, self.mismatch_count - 1)
                    self.target_change_count += 1
                return
        
    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()