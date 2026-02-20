"""
pcr/designers/qpcr.py
Real-time PCR (qPCR) 전용 Designer.
전략: Probe First -> Filter -> Primer Design around Probe (Config Driven)
"""
import primer3
from typing import Dict, Any, List, Optional

from .base import BasePrimerDesigner
from ..config.schema.root import BaseDesignInput, BaseDesignOutput
from ..components.primer import Primer, Probe
from ..components.amplicon import Amplicon

class QPCRPrimerDesigner(BasePrimerDesigner):
    """
    qPCR Designer (Probe-First Strategy)
    Config의 Tm Difference 설정을 기반으로 동적으로 Primer Tm 범위를 계산합니다.
    """
    ASSAY_TYPE = "TaqMan-qPCR"

    def __init__(self, input_data: BaseDesignInput):
        super().__init__(input_data)
        
        # Probe 주변 Gap 설정 (Config에 없다면 기본값 유지하되, 필요시 Config로 이동 가능)
        self.probe_gap = 3 

    def design(self) -> BaseDesignOutput:
        try:
            # Step 1: Probe 후보군 검색
            raw_probes = self._design_probes_only()

            if not raw_probes:
                return BaseDesignOutput(
                    amplicons=[], 
                    status="no_probes_found",
                    log_messages=["Primer3 found 0 probes."]
                )

            final_amplicons = []
            limit = self.input.overrides.get("PRIMER_NUM_RETURN", 100)
            processed_count = 0
            
            # Step 2: 각 Probe에 대해 Primer 매칭
            for probe_idx, probe_obj in enumerate(raw_probes):
                if processed_count >= limit:
                    break

                # 2-1. Filter
                if not self._is_valid_probe(probe_obj):
                    continue
                # 2-2. Primer Design (Config 기반 Tm 계산)
                amplicons = self._design_primers_for_probe(probe_obj, probe_idx)
                #print(amplicons)
                if amplicons:
                    final_amplicons.extend(amplicons)
                    processed_count += 1

            if not final_amplicons:
                return BaseDesignOutput(amplicons=[], status="no_valid_pairs")

            return BaseDesignOutput(
                amplicons=final_amplicons,
                status="success",
                log_messages=[
                    f"Raw Probes: {len(raw_probes)}",
                    f"Final Amplicons: {len(final_amplicons)}"
                ]
            )

        except Exception as e:
            import traceback
            traceback.print_exc()
            return BaseDesignOutput(amplicons=[], status="error", error_msg=str(e))

    def _design_probes_only(self) -> List[Probe]:
        """Probe Only Design"""
        probe_args = self.config.pcr_params.probe_kwargs.to_global_args()
        
        run_args = self.global_args.copy()
        run_args.update(probe_args)
        run_args.update({
            "PRIMER_PICK_LEFT_PRIMER": 0,
            "PRIMER_PICK_RIGHT_PRIMER": 0,
            "PRIMER_PICK_INTERNAL_OLIGO": 1,
            "PRIMER_INTERNAL_NUM_RETURN": 100,
        })

        res = primer3.bindings.design_primers(self.seq_args, run_args)
        
        num = res.get("PRIMER_INTERNAL_NUM_RETURNED", 0)
        probes = []
        for i in range(num):
            p = Probe.from_primer3(res, i)
            if p: probes.append(p)
        return probes

    def _is_valid_probe(self, probe: Probe) -> bool:
        """Design Constraints Filter"""
        seq = probe.sequence.upper()
        criteria = self.config.qc_criteria.probe
        p_start = probe.start_index
        p_len = len(probe.sequence)
        p_end = p_start + p_len # exclusive end

        t_start = self.input.target_start
        t_end = self.input.target_end

        # [조건] Probe 구간 안에 Target 구간이 완전히 들어와야 함
        # (Target이 Probe의 시작보다 뒤에 있고, 끝보다 앞에 있어야 함)
        if not (p_start <= t_start and p_end >= t_end):
            return False  # 칼같이 제거
        if criteria.avoid_5_prime_g and seq.startswith("G"):
            return False
        if "G" * (criteria.max_probe_poly_g + 1) in seq:
            return False
            
        return True

    def _design_primers_for_probe(self, probe: Probe, probe_idx: int) -> List[Amplicon]:
        """
        [수정됨] Config 값을 기반으로 Primer Tm 범위 동적 계산
        """
        # 1. Config에서 기준값 가져오기
        criteria = self.config.qc_criteria.probe
        min_diff = criteria.min_primer_probe_tm_diff
        max_diff = criteria.max_primer_probe_tm_diff
        
        primer_max_tm = probe.tm - min_diff
        primer_min_tm = probe.tm - max_diff
        # Opt Tm은 그 중간값으로 설정
        primer_opt_tm = (primer_max_tm + primer_min_tm) / 2.0
        # 3. Primer3 설정 업데이트
        run_args = self.global_args.copy()
        run_args.update({
            "PRIMER_OPT_TM": primer_opt_tm,
            "PRIMER_MIN_TM": primer_min_tm,
            "PRIMER_MAX_TM": primer_max_tm,
            

            # 기존 70-150bp -> 60-300bp로 확장하여 'unacceptable product size' 에러 방지
            "PRIMER_PRODUCT_SIZE_RANGE": [[60, 300]], 
            
            # (C) GC 함량 완화
            # 까다로운 서열(GC content failed)을 위해 범위를 넓힘
            "PRIMER_MIN_GC": 30.0,
            "PRIMER_MAX_GC": 75.0,
            
            # (D) 기본 필수 설정
            "PRIMER_PICK_INTERNAL_OLIGO": 1,
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_NUM_RETURN": 10,
        })

        # 4. Seq Args (Probe Fix & Gap)
        current_seq_args = self.seq_args.copy()
        current_seq_args["SEQUENCE_INTERNAL_OLIGO"] = probe.sequence
        
        p_start = probe.start_index if probe.start_index is not None else self.input.template_sequence.find(probe.sequence)
        
        if p_start >= 0:
            p_len = len(probe.sequence)
            excl_start = max(0, p_start - self.probe_gap)
            excl_end = min(len(self.input.template_sequence), p_start + p_len + self.probe_gap)
            current_seq_args["SEQUENCE_EXCLUDED_REGION"] = [[excl_start, excl_end - excl_start]]
        
        # 5. 실행 및 결과 변환
        res = primer3.bindings.design_primers(current_seq_args, run_args)
        num = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        amplicons = []
        for i in range(num):
            fwd = Primer.from_primer3(res, i, "LEFT")
            rev = Primer.from_primer3(res, i, "RIGHT")
            res_probe = Probe.from_primer3(res, i)
            
            if fwd and rev and res_probe:
                pair_penalty = float(res.get(f'PRIMER_PAIR_{i}_PENALTY', 0.0))
                amp_id = f"{self.input.name}_P{probe_idx}_{i}"
                amp = Amplicon(
                    id=amp_id,
                    forward=fwd,
                    reverse=rev,
                    probe=res_probe,
                    template_sequence=self.input.template_sequence,
                    reference_sequence=self.input.reference_sequence,
                    target_start_index=self.input.target_start,
                    target_end_index=self.input.target_end,
                    reference_id=self.input.reference_name,
                    pair_penalty=pair_penalty + probe.penalty
                )
                amp.is_qc_pass = True
                amplicons.append(amp)
                
        return amplicons