"""
pcr/designers/base/designer.py
모든 PCR designer의 공통 추상 클래스.
- config 로드 + Primer3 args 준비는 여기서
- design() 실행 흐름은 서브클래스에서 override
"""
import primer3
from typing import Dict, Any, List
from abc import ABC, abstractmethod

from .schema import BaseDesignInput, BaseDesignOutput
from pcr.components.primer import Primer, Probe
from pcr.components.amplicon import Amplicon


class BasePrimerDesigner(ABC):
    """
    공통 흐름:
      __init__  : config 로드 → assay별 overwrite(_apply_assay_config) → Primer3 args 준비
      design()  : Primer3 실행 → _process_results → 반환   (서브클래스에서 override 가능)
    """
    ASSAY_TYPE: str = "Base-PCR"

    def __init__(self, input_data: BaseDesignInput):
        self.input  = input_data
        self.config = input_data.config          # AppConfig (base값 이미 포함)
        self.design_name = f"{self.input.name}_{self.ASSAY_TYPE}"

        # ── assay별 config overwrite (서브클래스에서 구현) ──
        self._apply_assay_config()

        # ── Primer3 args 준비 ──
        self.seq_args: Dict[str, Any]    = {}
        self.global_args: Dict[str, Any] = {}
        self._prepare_primer3_args()
        
    @staticmethod
    def _reverse_complement(seq: str) -> str:
        return seq.translate(str.maketrans('ATGCatgcNn', 'TACGtacgNn'))[::-1]

    # ------------------------------------------------------------------
    # 서브클래스 구현 포인트
    # ------------------------------------------------------------------
    def _apply_assay_config(self):
        """
        assay별 파라미터를 self.config에 직접 overwrite.
        base는 no-op; 각 designer에서 필요한 값만 덮어씀.
        """
        pass

    @staticmethod
    def parse_sequence_with_brackets(raw_seq: str):
        """
        대괄호 표기법을 사용하여 순수 서열과 타겟 위치(1-based)를 반환합니다.
        예: 'ATGC[CG]ATGC' (MS-PCR) -> 순수 서열, [5, 6] 반환
        예: 'ATGC[A]TGC' (AS-PCR 변이) -> 순수 서열, [5] 반환
        """
        seq = raw_seq.replace(" ", "").replace("\n", "").upper()
        clean_seq = ""
        target_indices = []
        
        in_target = False
        clean_idx = 0
        
        for char in seq:
            if char == '[':
                in_target = True
            elif char == ']':
                in_target = False
            else:
                # 🔥 'C' 검사 조건을 없애고, 괄호 안의 '모든' 염기 위치를 타겟으로 등록!
                if in_target:
                    target_indices.append(clean_idx + 1) 
                
                clean_seq += char
                clean_idx += 1
                
        return clean_seq, target_indices
    
    @abstractmethod
    def design(self) -> BaseDesignOutput:
        """메인 파이프라인. 서브클래스에서 반드시 구현."""
        ...

    # ------------------------------------------------------------------
    # 공통 유틸
    # ------------------------------------------------------------------
    def _prepare_primer3_args(self):
        target_len = self.input.target_end - self.input.target_start

        self.seq_args = {
            'SEQUENCE_ID':       self.design_name,
            'SEQUENCE_TEMPLATE': self.input.template_sequence,
            'SEQUENCE_TARGET':   [self.input.target_start, target_len],
        }

        self.global_args = self.config.pcr_params.primer_kwargs.to_global_args()

        if self.config.pcr_params.probe_kwargs:
            # [MODIFIED] 변경 이유: 이전 논의대로, probe의 절대 Tm을 primer Tm 기준으로 동적 계산하기 위해 primer_kwargs를 인자로 전달
            probe_args = self.config.pcr_params.probe_kwargs.to_global_args(
                self.config.pcr_params.primer_kwargs
            )
            self.global_args.update(probe_args)
            self.global_args['PRIMER_PICK_INTERNAL_OLIGO'] = 1
        else:
            self.global_args['PRIMER_PICK_INTERNAL_OLIGO'] = 0

        if self.input.overrides:
            self.global_args.update(self.input.overrides)

    def _process_results(self, result: Dict[str, Any]) -> List[Amplicon]:
        """Primer3 raw 결과 → Amplicon 리스트"""
        num = result.get('PRIMER_PAIR_NUM_RETURNED', 0)
        amplicons = []

        for i in range(num):
            fwd = Primer.from_primer3(result, i, "LEFT")
            rev = Primer.from_primer3(result, i, "RIGHT")
            if not fwd or not rev:
                continue

            probe = None
            if self.config.pcr_params.probe_kwargs and \
               f"PRIMER_INTERNAL_{i}_SEQUENCE" in result:
                probe = Probe.from_primer3(result, i)

            amp = Amplicon(
                id=f"{self.input.name}_{i}",
                forward=fwd,
                reverse=rev,
                probe=probe,
                template_sequence=self.input.template_sequence,
                reference_sequence=self.input.reference_sequence or "",
                target_start_index=self.input.target_start,
                target_end_index=self.input.target_end,
                reference_id=self.input.reference_name,
                pair_penalty=float(result.get(f'PRIMER_PAIR_{i}_PENALTY', 0.0)),
            )
            # [MODIFIED] 임시로 True 할당 (이후 QCExecutor에서 판단)
            amp.is_qc_pass = True
            amplicons.append(amp)

        return amplicons

    def _run_primer3(self) -> Dict[str, Any]:
        return primer3.bindings.design_primers(self.seq_args, self.global_args)