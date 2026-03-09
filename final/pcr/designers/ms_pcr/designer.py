"""
pcr/designers/ms_pcr/designer.py
MS-PCR 전용 디자이너 클래스.
BasePrimerDesigner를 상속받아 M-Allele과 U-Allele에 대해 각각 프라이머를 설계하고 결과를 병합합니다.
"""
import copy
from typing import Dict, Any, List

import primer3

from ..base.designer import BasePrimerDesigner
from ..base.schema import BaseDesignInput, BaseDesignOutput
from pcr.components.primer import Primer
from pcr.components.amplicon import Amplicon
from pcr.components.region import GenomicRegion

class MSPCRPrimerDesigner(BasePrimerDesigner):
    ASSAY_TYPE = "MSPCR"

    def __init__(self, input_data: BaseDesignInput):
        super().__init__(input_data)
        
        # MS-PCR 전용 추가 입력값 로드
        self.templates = getattr(input_data, "templates", {})
        if "M" not in self.templates or "U" not in self.templates:
            raise ValueError("MS-PCR requires both 'M' and 'U' templates.")
            
        # Target CpG 좌표 (1-based 기준)
        self.target_cpgs = getattr(input_data, "target_cpg_indices", [])

    def _apply_assay_config(self):
        """MS-PCR은 프로브를 사용하지 않고, 프라이머 3' 말단에 집중합니다."""
        if hasattr(self.config.pcr_params, "probe_kwargs"):
            self.config.pcr_params.probe_kwargs = None

        # MS-PCR에서는 3' 말단이 중요하므로 3' 말단 GC 제약 등을 완화하거나 조절할 수 있습니다.
        # (필요 시 self.config.pcr_params.primer_kwargs의 속성을 변경)

    def _run_engine_for_allele(self, allele_type: str, seq: str) -> Dict[str, Any]:
        """특정 Allele (M 또는 U) 템플릿에 대해 Primer3를 실행합니다."""
        
        # 깊은 복사로 기존 인자 보호
        local_seq_args = copy.deepcopy(self.seq_args)
        local_global_args = copy.deepcopy(self.global_args)

        # 현재 알렐의 서열로 교체
        local_seq_args['SEQUENCE_TEMPLATE'] = seq
        
        # 내부 프로브 비활성화 강제
        local_global_args['PRIMER_PICK_INTERNAL_OLIGO'] = 0

        # [핵심] MS-PCR은 CpG 위치에 프라이머가 걸쳐야 함
        # 현재는 Target(CpG 포함 구간)을 감싸는(Flanking) 증폭을 수행함
        # 더 엄격한 MS-PCR의 경우 프라이머 3' 말단이 Target CpG에 위치하도록 강제할 수 있음
        
        try:
            return primer3.bindings.design_primers(local_seq_args, local_global_args)
        except Exception as e:
            return {"PRIMER_ERROR": str(e), "PRIMER_PAIR_NUM_RETURNED": 0}

    def _build_alignment_visual(self, m_amp: Amplicon, u_amp: Amplicon) -> List[str]:
        """M과 U 알렐의 프라이머 서열 정렬 상태를 텍스트로 시각화합니다."""
        lines = []
        
        raw_seq = self.input.template_sequence
        m_seq = self.templates.get("M", "")
        u_seq = self.templates.get("U", "")

        # 가장 넓은 증폭 구간 찾기 (시각화 범위)
        min_start = min(m_amp.forward.start_index, u_amp.forward.start_index)
        max_end = max(m_amp.reverse.end_index, u_amp.reverse.end_index)

        lines.append(f"RAW_SEQ  : {raw_seq[min_start:max_end]}")
        lines.append(f"M_ALLELE : {m_seq[min_start:max_end]}")
        lines.append(f"U_ALLELE : {u_seq[min_start:max_end]}")
        lines.append("-" * (max_end - min_start + 11))
        
        # 공백 패딩을 맞추기 위한 함수
        def pad(start_idx, seq):
            return " " * (start_idx - min_start) + seq

        lines.append(f"M_FWD    : {pad(m_amp.forward.start_index, m_amp.forward.sequence)}")
        lines.append(f"U_FWD    : {pad(u_amp.forward.start_index, u_amp.forward.sequence)}")
        
        m_rev_seq_rc = self._reverse_complement(m_amp.reverse.sequence)
        u_rev_seq_rc = self._reverse_complement(u_amp.reverse.sequence)
        
        # Reverse 프라이머는 3'->5' 방향이므로 RC 상태로 렌더링
        lines.append(f"M_REV_RC : {pad(m_amp.reverse.start_index, m_rev_seq_rc)}")
        lines.append(f"U_REV_RC : {pad(u_amp.reverse.start_index, u_rev_seq_rc)}")
        
        return lines

    def design(self) -> BaseDesignOutput:
        # 1. M-Allele에 대한 설계
        m_result = self._run_engine_for_allele("M", self.templates["M"])
        # 2. U-Allele에 대한 설계
        u_result = self._run_engine_for_allele("U", self.templates["U"])

        m_num = m_result.get('PRIMER_PAIR_NUM_RETURNED', 0)
        u_num = u_result.get('PRIMER_PAIR_NUM_RETURNED', 0)

        flat_amplicons: List[Amplicon] = []

        if m_num == 0 and u_num == 0:
            error_msg = f"Primer3 Failed for both M and U alleles. M_Explain: {m_result.get('PRIMER_LEFT_EXPLAIN', '')}"
            return BaseDesignOutput(status="fail", error_msg=error_msg, amplicons=[])

        # 3. M과 U 결과를 순위(Rank) 기준으로 짝짓기
        # MS-PCR에서는 M용과 U용 프라이머 쌍이 물리적으로 유사한 위치에서 작동하는 것이 이상적이므로
        # 같은 Rank(Primer3가 반환한 순서)끼리 묶어서 하나의 Set로 취급합니다.
        min_pairs = min(m_num, u_num)
        
        chrom = self.input.reference_name
        offset = getattr(self.input, "template_genomic_start", 0)

        for i in range(min_pairs):
            set_id = f"Set_{i}"
            
            # --- M Allele Amplicon 생성 ---
            m_fwd = Primer.from_primer3(m_result, i, "LEFT")
            m_rev = Primer.from_primer3(m_result, i, "RIGHT")
            
            if m_fwd and m_rev:
                m_fwd.template_sequence = self.templates["M"]
                m_rev.template_sequence = self.templates["M"]
                m_pair_penalty = float(m_result.get(f'PRIMER_PAIR_{i}_PENALTY', 0.0))
                
                # 절대 좌표 부여 (옵션)
                if offset > 0:
                    m_fwd.region = GenomicRegion(chrom, offset + m_fwd.start_index + 1, offset + m_fwd.end_index, "+")
                    m_rev.region = GenomicRegion(chrom, offset + m_rev.start_index + 1, offset + m_rev.end_index, "-")
                    
                m_amp = Amplicon(
                    id=f"{self.input.name}_{set_id}_M",
                    forward=m_fwd, reverse=m_rev, probe=None,
                    template_sequence=self.templates["M"],
                    target_start_index=self.input.target_start,
                    target_end_index=self.input.target_end,
                    pair_penalty=m_pair_penalty,
                )
                m_amp.allele_type = "M"
                m_amp.set_id = set_id
                m_amp.is_qc_pass = True

            # --- U Allele Amplicon 생성 ---
            u_fwd = Primer.from_primer3(u_result, i, "LEFT")
            u_rev = Primer.from_primer3(u_result, i, "RIGHT")
            
            if u_fwd and u_rev:
                u_fwd.template_sequence = self.templates["U"]
                u_rev.template_sequence = self.templates["U"]
                u_pair_penalty = float(u_result.get(f'PRIMER_PAIR_{i}_PENALTY', 0.0))
                
                # 절대 좌표 부여 (옵션)
                if offset > 0:
                    u_fwd.region = GenomicRegion(chrom, offset + u_fwd.start_index + 1, offset + u_fwd.end_index, "+")
                    u_rev.region = GenomicRegion(chrom, offset + u_rev.start_index + 1, offset + u_rev.end_index, "-")

                u_amp = Amplicon(
                    id=f"{self.input.name}_{set_id}_U",
                    forward=u_fwd, reverse=u_rev, probe=None,
                    template_sequence=self.templates["U"],
                    target_start_index=self.input.target_start,
                    target_end_index=self.input.target_end,
                    pair_penalty=u_pair_penalty,
                )
                u_amp.allele_type = "U"
                u_amp.set_id = set_id
                u_amp.is_qc_pass = True

            # 양쪽 모두 정상 생성되었을 때만 Set로 인정
            if m_fwd and m_rev and u_fwd and u_rev:
                # 시각화 데이터 생성 (M과 U를 묶어서 비교)
                alignment_view = self._build_alignment_visual(m_amp, u_amp)
                m_amp.alignment_visual = alignment_view
                u_amp.alignment_visual = alignment_view
                
                flat_amplicons.append(m_amp)
                flat_amplicons.append(u_amp)

        return BaseDesignOutput(status="success", amplicons=flat_amplicons, metadata={"assay": self.ASSAY_TYPE})