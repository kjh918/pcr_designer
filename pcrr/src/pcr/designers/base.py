import primer3
from typing import Any, Dict, List, Optional
from ..components import Amplicon, Primer
from ..config.schema.pcr import PrimerKwargs, ProbeKwargs

# ----------------------------------------------------------------------
# 1. BasePrimerDesigner (공통 부모)
# ----------------------------------------------------------------------
class BasePrimerDesigner:
    def __init__(
        self,
        template_sequence: str,
        target_start_index: int,
        target_end_index: int,
        config: PrimerKwargs,
        overrides: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.template_sequence = template_sequence
        self.target_start_index = int(target_start_index)
        self.target_end_index = int(target_end_index)
        
        self.cfg = config
        self.overrides = overrides or {}

        # Primer3 Arguments (딕셔너리)
        self.primer3_seq_args: Dict[str, Any] = {}
        self.primer3_global_args: Dict[str, Any] = {}
        self.primer3_result: Optional[Dict[str, Any]] = None
        self.amplicon_list: List[Amplicon] = []

        # 초기 설정 로드
        self._init_primer3_args()

    def _init_primer3_args(self) -> None:
        target_len = self.target_end_index - self.target_start_index + 1
        self.primer3_seq_args = {
            "SEQUENCE_ID": "generic_design",
            "SEQUENCE_TEMPLATE": self.template_sequence,
            "SEQUENCE_TARGET": [self.target_start_index, target_len],
        }

        # Config 객체 -> 딕셔너리 변환
        if hasattr(self.cfg, "to_global_args"):
            self.primer3_global_args = self.cfg.to_global_args()
        else:
            self.primer3_global_args = self.cfg.to_primer3_args()

        # Overrides 적용
        if self.overrides:
            self.primer3_global_args.update(self.overrides)

    def run_primer3(self) -> None:
        """현재 설정된 args로 Primer3 실행"""
        try:
            self.primer3_result = primer3.bindings.design_primers(
                seq_args=self.primer3_seq_args,
                global_args=self.primer3_global_args,
            )
        except Exception:
            self.primer3_result = {}

    def reset(self) -> None:
        """설정 초기화 (재사용 시)"""
        self._init_primer3_args()
        self.primer3_result = None
        self.amplicon_list = []

    def design(self) -> List[Amplicon]:
        self.run_primer3()
        return self._build_amplicons()

    def _build_amplicons(self) -> List[Amplicon]:
        raise NotImplementedError


# ----------------------------------------------------------------------
# 2. PrimerDesigner (일반)
# ----------------------------------------------------------------------
class PrimerDesigner(BasePrimerDesigner):
    def _build_amplicons(self) -> List[Amplicon]:
        if not self.primer3_result: return []
        n_pairs = int(self.primer3_result.get("PRIMER_PAIR_NUM_RETURNED", 0))
        amplicons = []
        for i in range(n_pairs):
            # ... (기존 코드와 동일: Fwd/Rev Primer 및 Amplicon 생성) ...
            left_seq = self.primer3_result.get(f"PRIMER_LEFT_{i}_SEQUENCE")
            right_seq = self.primer3_result.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
            if not (left_seq and right_seq): continue
            
            fwd = Primer(
                template_sequence=self.template_sequence, sequence=left_seq,
                tm=self.primer3_result.get(f"PRIMER_LEFT_{i}_TM"), strand="forward", primer_type="forward_primer",
                target_start_index=self.target_start_index, target_end_index=self.target_end_index,
                binding_start_index=self.primer3_result.get(f"PRIMER_LEFT_{i}")[0],
                binding_end_index=self.primer3_result.get(f"PRIMER_LEFT_{i}")[0] + self.primer3_result.get(f"PRIMER_LEFT_{i}")[1]
            )
            rev_data = self.primer3_result.get(f"PRIMER_RIGHT_{i}")
            rev = Primer(
                template_sequence=self.template_sequence, sequence=right_seq,
                tm=self.primer3_result.get(f"PRIMER_RIGHT_{i}_TM"), strand="reverse", primer_type="reverse_primer",
                target_start_index=self.target_start_index, target_end_index=self.target_end_index,
                binding_start_index=rev_data[0] - rev_data[1] + 1, binding_end_index=rev_data[0] + 1
            )
            amplicons.append(Amplicon(
                template_sequence=self.template_sequence, target_start_index=self.target_start_index, target_end_index=self.target_end_index,
                forward_primer=fwd, reverse_primer=rev,
                product_size=self.primer3_result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
                tm=self.primer3_result.get(f"PRIMER_PAIR_{i}_PRODUCT_TM")
            ))
        return amplicons


# ----------------------------------------------------------------------
# 3. ProbePrimerDesigner (과거 코드의 유연성 + 현재 구조 통합)
# ----------------------------------------------------------------------
class ProbePrimerDesigner(BasePrimerDesigner):
    def __init__(
        self,
        template_sequence: str,
        target_start_index: int,
        target_end_index: int,
        config: PrimerKwargs,
        probe_config: Optional[ProbeKwargs] = None, # ✅ 추가됨
        overrides: Optional[Dict[str, Any]] = None,
    ) -> None:
        # ✅ 부모에게는 probe_config를 넘기지 않음 (TypeError 해결)
        super().__init__(
            template_sequence=template_sequence,
            target_start_index=target_start_index,
            target_end_index=target_end_index,
            config=config,
            overrides=overrides
        )
        self.probe_cfg = probe_config

    def _configure_probe_only(self) -> None:
        """Probe만 뽑도록 설정"""
        self.primer3_global_args["PRIMER_PICK_LEFT_PRIMER"] = 0
        self.primer3_global_args["PRIMER_PICK_RIGHT_PRIMER"] = 0
        self.primer3_global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 1
        
        # Probe Config 적용
        if self.probe_cfg:
            self.primer3_global_args.update(self.probe_cfg.to_global_args())
            # Target 강제 Overlap 로직 (제공해주신 코드)
            probe_min_len = self.probe_cfg.min_length
            
            excl_left_end = max(0, self.target_end_index - probe_min_len)
            excl_right_start = min(len(self.template_sequence), self.target_start_index + probe_min_len)
            
            exclusions = []
            if excl_left_end > 0: exclusions.append([0, excl_left_end])
            if len(self.template_sequence) - excl_right_start > 0:
                exclusions.append([excl_right_start, len(self.template_sequence) - excl_right_start])
            
            if exclusions:
                self.primer3_seq_args["SEQUENCE_INTERNAL_EXCLUDED_REGION"] = exclusions

        if self.overrides:
            self.primer3_global_args.update(self.overrides)

    def _design_primers_for_probe(self, probe_amp: Amplicon) -> List[Amplicon]:
        """고정된 Probe에 맞춰 Primer Pair 디자인 (Retry 로직 포함)"""
        probe = probe_amp.probe
        if not probe or not probe.tm: return []

        # 1. 초기화 및 Probe 고정
        self.reset()
        self.primer3_global_args["PRIMER_PICK_LEFT_PRIMER"] = 1
        self.primer3_global_args["PRIMER_PICK_RIGHT_PRIMER"] = 1
        self.primer3_global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 1
        self.primer3_seq_args["SEQUENCE_INTERNAL_OLIGO"] = probe.sequence
        
        # Probe 디자인 시 사용했던 Target 제약 제거 (Primer는 자유롭게)
        if "SEQUENCE_TARGET" in self.primer3_seq_args:
            del self.primer3_seq_args["SEQUENCE_TARGET"]
        
        # Product Size Range 확장 (필수)
        #self.primer3_global_args["PRIMER_PRODUCT_SIZE_RANGE"] = [[80, len(self.template_sequence) - 10]]

        # 2. 설정값 준비
        gap = self.overrides.get("min_primer_probe_distance", 5)
        target_tm = probe.tm
        
        if self.probe_cfg:
            min_diff = self.probe_cfg.min_primer_probe_tm_diff
            max_diff = self.probe_cfg.max_primer_probe_tm_diff
        else:
            min_diff = 5.0
            max_diff = 10.0
        
        min_diff = self.overrides.get("min_primer_probe_tm_diff", min_diff)
        max_diff = self.overrides.get("max_primer_probe_tm_diff", max_diff)

        # 3. 1차 시도 (Strict)
        results = self._run_primer3_with_params(probe, gap, target_tm, min_diff, max_diff)
        return results

        ## 4. 2차 시도 (Relaxed) - 실패 시 조건 완화
        #relaxed_gap = 1
        #relaxed_min_diff = 2.0
        #relaxed_max_diff = 15.0
        #return self._run_primer3_with_params(probe, relaxed_gap, target_tm, relaxed_min_diff, relaxed_max_diff)

    def _run_primer3_with_params(self, probe, gap, target_tm, min_diff, max_diff) -> List[Amplicon]:
        # Exclusion Zone
        p_start = probe.binding_start_index
        p_len = len(probe.sequence)
        excl_start = max(0, p_start - gap)
        excl_end = min(len(self.template_sequence), p_start + p_len + gap)
        excl_len = excl_end - excl_start
        if excl_len > 0:
            self.primer3_seq_args["SEQUENCE_EXCLUDED_REGION"] = [[excl_start, excl_len]]
        
        # Tm Window
        self.primer3_global_args["PRIMER_OPT_TM"] = (target_tm - min_diff + target_tm - max_diff) / 2
        self.primer3_global_args["PRIMER_MAX_TM"] = target_tm - min_diff
        self.primer3_global_args["PRIMER_MIN_TM"] = target_tm - max_diff
        
        self.run_primer3()
        
        if not self.primer3_result: return []
        
        # 파싱
        n_pairs = int(self.primer3_result.get("PRIMER_PAIR_NUM_RETURNED", 0))
        results = []
        for i in range(n_pairs):
            left_seq = self.primer3_result.get(f"PRIMER_LEFT_{i}_SEQUENCE")
            right_seq = self.primer3_result.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
            left_loc = self.primer3_result.get(f"PRIMER_LEFT_{i}")
            right_loc = self.primer3_result.get(f"PRIMER_RIGHT_{i}")
            
            if not (left_seq and right_seq and left_loc and right_loc): continue

            fwd = Primer(
                template_sequence=self.template_sequence, sequence=left_seq,
                tm=self.primer3_result.get(f"PRIMER_LEFT_{i}_TM"), strand="forward", primer_type="forward_primer",
                target_start_index=self.target_start_index, target_end_index=self.target_end_index,
                binding_start_index=left_loc[0], binding_end_index=left_loc[0] + left_loc[1]
            )
            r_start, r_len = right_loc
            rev = Primer(
                template_sequence=self.template_sequence, sequence=right_seq,
                tm=self.primer3_result.get(f"PRIMER_RIGHT_{i}_TM"), strand="reverse", primer_type="reverse_primer",
                target_start_index=self.target_start_index, target_end_index=self.target_end_index,
                binding_start_index=r_start - r_len + 1, binding_end_index=r_start + 1
            )
            results.append(Amplicon(
                template_sequence=self.template_sequence, target_start_index=self.target_start_index, target_end_index=self.target_end_index,
                forward_primer=fwd, reverse_primer=rev, probe=probe,
                product_size=self.primer3_result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
                tm=self.primer3_result.get(f"PRIMER_PAIR_{i}_PRODUCT_TM")
            ))
        return results

    def design(self) -> List[Amplicon]:
        # 1. Probe Only
        self.reset()
        self._configure_probe_only()
        self.run_primer3()
        
        # Probe 파싱 (필터 적용)
        probe_candidates = []
        if self.primer3_result:
            n_probes = int(self.primer3_result.get("PRIMER_INTERNAL_NUM_RETURNED", 0))
            
            # 필터 설정
            if self.probe_cfg:
                max_poly_g = self.probe_cfg.max_probe_poly_g
                max_3prime_gc = self.probe_cfg.max_probe_3_end_gc
            else:
                max_poly_g = 4; max_3prime_gc = 2

            for i in range(n_probes):
                seq = self.primer3_result.get(f"PRIMER_INTERNAL_{i}_SEQUENCE")
                if not seq: continue
                # 필터링
                if seq.startswith("G"): continue
                if "G" * (max_poly_g + 1) in seq: continue
                tail_5bp = seq[-5:]
                if tail_5bp.count("G") + tail_5bp.count("C") >= max_3prime_gc: continue
                
                tm = self.primer3_result.get(f"PRIMER_INTERNAL_{i}_TM")
                probe = Primer(
                    template_sequence=self.template_sequence, sequence=seq,
                    tm=tm, strand="forward", primer_type="probe",
                    target_start_index=self.target_start_index, target_end_index=self.target_end_index,
                )
                probe_candidates.append(Amplicon(
                    template_sequence=self.template_sequence, target_start_index=self.target_start_index, target_end_index=self.target_end_index,
                    forward_primer=None, reverse_primer=None, probe=probe
                ))

        # 2. Primer Design loop
        final_results = []
        for probe_amp in probe_candidates:
            paired = self._design_primers_for_probe(probe_amp)
            final_results.extend(paired)
            
        self.amplicon_list = final_results
        return final_results

    def _build_amplicons(self) -> List[Amplicon]:
        return []