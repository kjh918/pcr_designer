"""
pcr/designers/ms_pcr/designer.py
MS-PCR 전용 디자이너 클래스.
M-Allele은 타겟 CpG에 Forward 프라이머의 3' 말단을 정확히 고정하고, 
U-Allele은 타겟 CpG가 Forward 프라이머의 3' 말단(default 3bp 이내)에 위치하도록 강제한 후, 
가장 물리적 위치가 일치하고 열역학적 점수가 좋은 M과 U를 묶어 세트로 구성합니다.
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
        
        self.templates = getattr(input_data, "templates", {})
        if "M" not in self.templates or "U" not in self.templates:
            raise ValueError("MS-PCR requires both 'M' and 'U' templates.")
            
        # Target CpG 좌표 (1-based 기준)
        self.target_cpgs = getattr(input_data, "target_cpg_indices", [])

    def _apply_assay_config(self):
        """MS-PCR은 프로브를 사용하지 않고, 프라이머 3' 말단에 집중합니다."""
        if hasattr(self.config.pcr_params, "probe_kwargs"):
            self.config.pcr_params.probe_kwargs = None

    def _run_engine_for_allele(self, allele_type: str, seq: str, cpg_idx: int, window_size: int = 1) -> List[Dict[str, Any]]:
        """
        Allele(M 또는 U)에 대해 프라이머를 탐색합니다.
        조건: Forward 프라이머의 3' 말단에서 최대 `window_size` 염기 이내에 타겟 CpG가 존재해야 합니다.
        - M의 경우 window_size=1로 고정하여 말단을 정확히 일치시킵니다.
        - U의 경우 window_size=3 등으로 여유를 주어 말단 근처에 타겟이 포함되도록 합니다.
        """
        local_seq_args = copy.deepcopy(self.seq_args)
        local_global_args = copy.deepcopy(self.global_args)
        
        local_seq_args['SEQUENCE_TEMPLATE'] = seq
        for key in ["SEQUENCE_TARGET", "SEQUENCE_FORCE_LEFT_START", "SEQUENCE_FORCE_RIGHT_START", "SEQUENCE_FORCE_RIGHT_END", "SEQUENCE_EXCLUDED_REGION"]:
            local_seq_args.pop(key, None)

        primer_kw = self.config.pcr_params.primer_kwargs
        n_primers = int(getattr(primer_kw, "n_candidates", 30))
        
        local_global_args.update({
            "PRIMER_PICK_LEFT_PRIMER": 1, 
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_PRIMER": 0,
            "PRIMER_NUM_RETURN": n_primers,
            "PRIMER_EXPLAIN_FLAG": 1,
        })
        
        candidates = []
        
        # 3' 말단을 타겟 위치(cpg_idx)부터 cpg_idx + window_size - 1 까지 이동하며 탐색
        for shift in range(window_size):
            end_pos = cpg_idx + shift
            if end_pos >= len(seq): break
            
            # Forward 프라이머 3' 말단 위치 강제 고정
            local_seq_args['SEQUENCE_FORCE_LEFT_END'] = end_pos
            
            try:
                res = primer3.bindings.design_primers(local_seq_args, local_global_args)
                print(res)
                num_ret = res.get('PRIMER_PAIR_NUM_RETURNED', 0)
                for i in range(num_ret):
                    fwd = Primer.from_primer3(res, i, "LEFT")
                    rev = Primer.from_primer3(res, i, "RIGHT")
                    
                    if fwd and rev:
                        penalty = float(res.get(f'PRIMER_PAIR_{i}_PENALTY', 9999.0))
                        uid = f"{fwd.start_index}_{rev.start_index}"
                        
                        candidates.append({
                            "uid": uid,
                            "fwd": fwd,
                            "rev": rev,
                            "penalty": penalty,
                            "shift": shift # 타겟이 말단에서 몇 칸 떨어져 있는지 기록
                        })
            except Exception:
                continue

        # Primer3 점수(Penalty) 오름차순으로 1차 정렬
        candidates.sort(key=lambda x: x["penalty"])
        
        # 중복된 프라이머 쌍(위치 기준) 제거
        unique_candidates = []
        seen = set()
        for c in candidates:
            if c["uid"] not in seen:
                seen.add(c["uid"])
                unique_candidates.append(c)

        # 상위 n_primers 개만 반환
        return unique_candidates[:n_primers]

    def _build_alignment_visual(self, m_amp: Amplicon, u_amp: Amplicon) -> List[str]:
        """M과 U 알렐의 프라이머 서열 정렬 상태를 텍스트로 시각화합니다."""
        lines = []
        raw_seq = self.input.template_sequence
        m_seq = self.templates.get("M", "")
        u_seq = self.templates.get("U", "")

        min_start = min(m_amp.forward.start_index, u_amp.forward.start_index)
        max_end = max(m_amp.reverse.end_index, u_amp.reverse.end_index)

        lines.append(f"RAW_SEQ  : {raw_seq[min_start:max_end]}")
        lines.append(f"M_ALLELE : {m_seq[min_start:max_end]}")
        lines.append(f"U_ALLELE : {u_seq[min_start:max_end]}")
        lines.append("-" * (max_end - min_start + 11))
        
        def pad(start_idx, seq):
            return " " * (start_idx - min_start) + seq

        lines.append(f"M_FWD    : {pad(m_amp.forward.start_index, m_amp.forward.sequence)}")
        lines.append(f"U_FWD    : {pad(u_amp.forward.start_index, u_amp.forward.sequence)}")
        
        m_rev_seq_rc = self._reverse_complement(m_amp.reverse.sequence)
        u_rev_seq_rc = self._reverse_complement(u_amp.reverse.sequence)
        
        lines.append(f"M_REV_RC : {pad(m_amp.reverse.start_index, m_rev_seq_rc)}")
        lines.append(f"U_REV_RC : {pad(u_amp.reverse.start_index, u_rev_seq_rc)}")
        
        return lines

    def design(self) -> BaseDesignOutput:
        if not self.target_cpgs:
            return BaseDesignOutput(status="fail", error_msg="Target CpG [CG] is missing.", amplicons=[])

        # 0-based 인덱스 (C의 위치)
        cpg_idx = self.target_cpgs[0] - 1
        chrom = self.input.reference_name
        offset = getattr(self.input, "template_genomic_start", 0)

        # U-Allele을 위한 3' 말단 윈도우 사이즈 설정 (API 입력이 없으면 기본 3으로 설정)
        u_window_size = getattr(self.input, "window_size_3prime", 3)

        # Step 1: M과 U에 대해 독립적으로 최적의 후보군 탐색
        # 🔥 M은 윈도우를 1로 고정하여 3' 말단에 정확히 일치시킴
        m_candidates = self._run_engine_for_allele("M", self.templates["M"], cpg_idx, window_size=1)
        # 🔥 U는 지정된 윈도우 사이즈(기본 3)만큼 탐색
        u_candidates = self._run_engine_for_allele("U", self.templates["U"], cpg_idx, window_size=u_window_size)
        print(u_candidates)
        if not m_candidates or not u_candidates:
            error_msg = f"Failed to design primers. Found M: {len(m_candidates)}, U: {len(u_candidates)}. Try relaxing GC/Tm constraints."
            return BaseDesignOutput(status="fail", error_msg=error_msg, amplicons=[])

        # Step 2: Cross-Pairing (M과 U 조합 묶기)
        cross_pairs = []
        for m in m_candidates:
            for u in u_candidates:
                # M과 U 프라이머의 위치(Index) 차이 계산
                fwd_diff = abs(m["fwd"].start_index - u["fwd"].start_index)
                rev_diff = abs(m["rev"].start_index - u["rev"].start_index)
                
                # [점수 계산식] 열역학적 페널티 합산 + 프라이머 결합 위치 차이에 대한 페널티 가중치
                # 위치가 비슷할수록(증폭 산물 크기가 같을수록) 좋은 세트로 평가됨
                combined_score = m["penalty"] + u["penalty"] + (fwd_diff * 0.5) + (rev_diff * 0.5)
                
                cross_pairs.append({
                    "m": m,
                    "u": u,
                    "score": combined_score
                })

        # 점수가 가장 좋은 순(오름차순)으로 정렬
        cross_pairs.sort(key=lambda x: x["score"])

        # Step 3: 완벽한 세트 조립
        flat_amplicons: List[Amplicon] = []
        primer_kw = self.config.pcr_params.primer_kwargs
        max_sets = int(getattr(primer_kw, "n_candidates", 50))
        
        seen_combinations = set()
        set_count = 0

        for pair in cross_pairs:
            if set_count >= max_sets:
                break
                
            m = pair["m"]
            u = pair["u"]
            
            # 이미 조립된 적 있는 프라이머 위치 조합은 건너뜀
            uid = f"M{m['fwd'].start_index}_{m['rev'].start_index}_U{u['fwd'].start_index}_{u['rev'].start_index}"
            if uid in seen_combinations:
                continue
            seen_combinations.add(uid)
            
            set_id = f"Set_{set_count}"
            set_count += 1
            
            # --- M Allele Amplicon 조립 ---
            m_fwd, m_rev = copy.deepcopy(m["fwd"]), copy.deepcopy(m["rev"])
            m_fwd.template_sequence, m_rev.template_sequence = self.templates["M"], self.templates["M"]
            
            if offset > 0:
                m_fwd.region = GenomicRegion(chrom, offset + m_fwd.start_index + 1, offset + m_fwd.end_index, "+")
                m_rev.region = GenomicRegion(chrom, offset + m_rev.start_index + 1, offset + m_rev.end_index, "-")
                
            m_amp = Amplicon(
                id=f"{self.input.name}_{set_id}_M",
                set_id=f"{self.input.name}_{set_id}", 
                forward=m_fwd, reverse=m_rev, probe=None,
                template_sequence=self.templates["M"],
                target_start_index=self.input.target_start, target_end_index=self.input.target_end,
                pair_penalty=m["penalty"],
            )
            m_amp.allele_type = "M"
            m_amp.is_qc_pass = True

            # --- U Allele Amplicon 조립 ---
            u_fwd, u_rev = copy.deepcopy(u["fwd"]), copy.deepcopy(u["rev"])
            u_fwd.template_sequence, u_rev.template_sequence = self.templates["U"], self.templates["U"]
            
            if offset > 0:
                u_fwd.region = GenomicRegion(chrom, offset + u_fwd.start_index + 1, offset + u_fwd.end_index, "+")
                u_rev.region = GenomicRegion(chrom, offset + u_rev.start_index + 1, offset + u_rev.end_index, "-")

            u_amp = Amplicon(
                id=f"{self.input.name}_{set_id}_U",
                set_id=f"{self.input.name}_{set_id}",
                forward=u_fwd, reverse=u_rev, probe=None,
                template_sequence=self.templates["U"],
                target_start_index=self.input.target_start, target_end_index=self.input.target_end,
                pair_penalty=u["penalty"],
            )
            u_amp.allele_type = "U"
            u_amp.is_qc_pass = True

            alignment_view = self._build_alignment_visual(m_amp, u_amp)
            m_amp.alignment_visual, u_amp.alignment_visual = alignment_view, alignment_view
            
            flat_amplicons.append(m_amp)
            flat_amplicons.append(u_amp)

        return BaseDesignOutput(status="success", amplicons=flat_amplicons, metadata={"assay": self.ASSAY_TYPE})