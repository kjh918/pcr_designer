import primer3
from typing import Dict, Any, List

from ..base.designer import BasePrimerDesigner
from ..base.schema import BaseDesignInput, BaseDesignOutput
from pcr.components.primer import Primer, Probe
from pcr.components.amplicon import Amplicon

class QPCRPrimerDesigner(BasePrimerDesigner):
    ASSAY_TYPE = "qPCR"

    def __init__(self, input_data: BaseDesignInput):
        super().__init__(input_data)
        self.probe_gap = self.config.pcr_params.probe_kwargs.probe_gap

    # ------------------------------------------------------------------
    # [Override] qPCR 전용 config overwrite
    # ------------------------------------------------------------------
    def _apply_assay_config(self):
        if not self.config.pcr_params.probe_kwargs:
            raise ValueError("qPCR 설계에는 probe_kwargs 설정이 필수입니다.")

    # ------------------------------------------------------------------
    # [Override] 메인 파이프라인
    # ------------------------------------------------------------------
    def design(self) -> BaseDesignOutput:
        try:
            raw_probes = self._design_probes_only()
            if not raw_probes:
                return BaseDesignOutput(status="no_probes_found", log_messages=["Primer3 found 0 probes."])

            limit = self.input.overrides.get("return_top_k", self.config.pcr_params.primer_kwargs.n_candidates)
            final_amplicons = []
            processed = 0

            for probe_idx, probe_obj in enumerate(raw_probes):
                if processed >= limit:
                    break
                if not self._is_valid_probe(probe_obj):
                    continue
                
                amps = self._design_primers_for_probe(probe_obj, probe_idx)

                if amps:
                    final_amplicons.extend(amps)
                    processed += 1

            if not final_amplicons:
                return BaseDesignOutput(status="no_valid_pairs")

            return BaseDesignOutput(
                amplicons=final_amplicons,
                status="success",
                log_messages=[f"Raw Probes: {len(raw_probes)}", f"Final Amplicons: {len(final_amplicons)}"]
            )

        except Exception as e:
            import traceback; traceback.print_exc()
            return BaseDesignOutput(status="error", error_msg=str(e))

    # ------------------------------------------------------------------
    # 내부 헬퍼
    # ------------------------------------------------------------------
    def _design_probes_only(self) -> List[Probe]:
        # 💡 [MODIFIED] 스키마 변경점 반영: Primer의 opt_tm을 넘겨서 Probe Tm 범위를 동적으로 계산합니다!
        primer_opt_tm = self.config.pcr_params.primer_kwargs.opt_tm
        probe_args = self.config.pcr_params.probe_kwargs.to_global_args(primer_opt_tm=primer_opt_tm)
        
        run_args = self.global_args.copy()
        run_args.update(probe_args)
        
        n_candidates = self.config.pcr_params.probe_kwargs.n_candidates

        run_args.update({
            "PRIMER_PICK_LEFT_PRIMER":     0,
            "PRIMER_PICK_RIGHT_PRIMER":    0,
            "PRIMER_PICK_INTERNAL_OLIGO":  1,
            "PRIMER_INTERNAL_NUM_RETURN":  n_candidates,
        })
        res  = primer3.bindings.design_primers(self.seq_args, run_args)
        num  = res.get("PRIMER_INTERNAL_NUM_RETURNED", 0)
        return [p for i in range(num) if (p := Probe.from_primer3(res, i))]

    def _is_valid_probe(self, probe: Probe) -> bool:
        seq = probe.sequence.upper()
        criteria = getattr(self.config.qc_criteria, "probe", None)
        if not criteria: return True
        print(seq)
        p_start, p_end = probe.start_index, probe.start_index + len(seq)
        if not (p_start <= self.input.target_start and p_end >= self.input.target_end):
            return False
        print(1)
        if criteria.avoid_5_prime_g and seq.startswith("G"):
            return False
        if "G" * (criteria.max_probe_poly_g + 1) in seq:
            return False
        print('\n'+seq)
        return True
    
    def _design_primers_for_probe(self, probe: Probe, probe_idx: int) -> List[Amplicon]:
        primer_kw = self.config.pcr_params.primer_kwargs
        probe_kw = self.config.pcr_params.probe_kwargs
        
        min_diff = probe_kw.min_tm_diff
        max_diff = probe_kw.max_tm_diff
        
        # 실제 뽑힌 Probe의 Tm을 기준으로 Primer의 상/하한선을 역추적합니다.
        primer_max_tm = probe.tm - min_diff
        primer_min_tm = probe.tm - max_diff
        primer_opt_tm = (primer_max_tm + primer_min_tm) / 2.0

        # 🔍 [DEBUG 1] Probe Tm과 계산된 Primer Tm 범위 출력
        print(f"\n🔍 [DEBUG - Probe {probe_idx}] Probe Seq: {probe.sequence} (Tm: {probe.tm:.1f}°C)")
        print(f"   => 요구되는 Primer Tm 범위: {primer_min_tm:.1f}°C ~ {primer_max_tm:.1f}°C (Opt: {primer_opt_tm:.1f}°C)")

        run_args = self.global_args.copy()
        run_args.update({
            "PRIMER_OPT_TM":              primer_opt_tm,
            "PRIMER_MIN_TM":              primer_min_tm,
            "PRIMER_MAX_TM":              primer_max_tm,
            "PRIMER_PRODUCT_SIZE_RANGE":  [[primer_kw.min_amplicon_length, primer_kw.max_amplicon_length]],
            "PRIMER_MIN_GC":              primer_kw.min_gc,
            "PRIMER_MAX_GC":              primer_kw.max_gc,
            "PRIMER_PICK_INTERNAL_OLIGO": 1,
            "PRIMER_PICK_LEFT_PRIMER":    1,
            "PRIMER_PICK_RIGHT_PRIMER":   1,
            "PRIMER_NUM_RETURN":          probe_kw.n_primers_per_probe,
        })

        seq_args = self.seq_args.copy()
        seq_args["SEQUENCE_INTERNAL_OLIGO"] = probe.sequence

        p_start = probe.start_index if probe.start_index is not None \
                  else self.input.template_sequence.find(probe.sequence)
        
        if p_start >= 0:
            excl_start = max(0, p_start - self.probe_gap)
            excl_end   = min(len(self.input.template_sequence), p_start + len(probe.sequence) + self.probe_gap)
            seq_args["SEQUENCE_EXCLUDED_REGION"] = [[excl_start, excl_end - excl_start]]
            
            # 🔍 [DEBUG 2] Excluded Region (프라이머가 오면 안 되는 구역) 출력
            print(f"   => Excluded Region (Probe위치): start={excl_start}, length={excl_end - excl_start}")

        # Primer3 실행
        res = primer3.bindings.design_primers(seq_args, run_args)
        num = res.get("PRIMER_PAIR_NUM_RETURNED", 0)

        # 🔍 [DEBUG 3] 실패 원인 분석 (EXPLAIN 출력)
        if num == 0:
            print(f"❌ [DEBUG - Probe {probe_idx}] Primer3가 프라이머 쌍을 찾지 못했습니다!")
            print(f"   [LEFT_EXPLAIN]  : {res.get('PRIMER_LEFT_EXPLAIN', '없음')}")
            print(f"   [RIGHT_EXPLAIN] : {res.get('PRIMER_RIGHT_EXPLAIN', '없음')}")
            print(f"   [PAIR_EXPLAIN]  : {res.get('PRIMER_PAIR_EXPLAIN', '없음')}")
        else:
            print(f"✅ [DEBUG - Probe {probe_idx}] 프라이머 {num}쌍 생성 성공!")

        amplicons = []
        for i in range(num):
            fwd, rev, prb = Primer.from_primer3(res, i, "LEFT"), Primer.from_primer3(res, i, "RIGHT"), Probe.from_primer3(res, i)
            if fwd and rev and prb:
                amp = Amplicon(
                    id=f"{self.input.name}_P{probe_idx}_{i}",
                    forward=fwd, reverse=rev, probe=prb,
                    template_sequence=self.input.template_sequence,
                    reference_sequence=self.input.reference_sequence or "",
                    target_start_index=self.input.target_start, target_end_index=self.input.target_end,
                    reference_id=self.input.reference_name,
                    pair_penalty=float(res.get(f'PRIMER_PAIR_{i}_PENALTY', 0.0)) + probe.penalty,
                )
                amp.is_qc_pass = True
                amplicons.append(amp)

        return amplicons