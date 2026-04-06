import primer3
from typing import Dict, Any, List, Tuple
import copy

from pcr.designers.base.designer import BasePrimerDesigner
from pcr.designers.base.schema import BaseDesignInput, BaseDesignOutput
from pcr.components.primer import Primer
from pcr.components.amplicon import Amplicon
from pcr.components.region import GenomicRegion

class ASPCRPrimerDesigner(BasePrimerDesigner):
    ASSAY_TYPE = "ASPCR"

    def __init__(self, input_data: BaseDesignInput):
        super().__init__(input_data)
        self.templates = getattr(input_data, "templates", {})
        
        # 🔥 HTML UI의 "+", "-" 값을 "forward", "reverse"로 자동 변환
        raw_fixed_prime = str(getattr(input_data, "fixed_prime", "forward")).strip().lower()
        if raw_fixed_prime == "+":
            self.fixed_prime = "forward"
        elif raw_fixed_prime == "-":
            self.fixed_prime = "reverse"
        else:
            self.fixed_prime = raw_fixed_prime

        self.target_index = input_data.target_start # SNP 위치 (0-based 기준점)
        
        # 🔥 HTML UI의 "None" 문자열을 안전하게 파싱 (None으로 처리)
        raw_mm_pos = getattr(input_data, "mismatch_pos", 3)
        if str(raw_mm_pos).strip().lower() == "none" or not raw_mm_pos:
            self.mismatch_pos = None
        else:
            self.mismatch_pos = int(raw_mm_pos)
            
        self.mismatch_intensity = getattr(input_data, "mismatch_intensity", "strong")

    # =========================================================================
    # Helpers (Coordinate-safe, Template-slice based)
    # =========================================================================
    @staticmethod
    def _rc(seq: str) -> str:
        return seq.translate(str.maketrans('ATGCatgcNn', 'TACGtacgNn'))[::-1]

    @staticmethod
    def _left_pos_to_span(pos: Tuple[int, int]) -> Tuple[int, int]:
        """Primer3 LEFT position (start, length) -> (start, end inclusive)"""
        start, length = int(pos[0]), int(pos[1])
        return start, start + length - 1

    @staticmethod
    def _right_pos_to_span(pos: Tuple[int, int]) -> Tuple[int, int]:
        """Primer3 RIGHT position (three_prime, length) -> (start, end inclusive)"""
        three_prime, length = int(pos[0]), int(pos[1])
        return three_prime - length + 1, three_prime

    @staticmethod
    def _validate_anchor(fixed_prime: str, left_span: Tuple[int, int], right_span: Tuple[int, int], target_index: int) -> bool:
        """Primer3가 target_index에 3' 말단을 정확히 고정했는지 재검증합니다."""
        if fixed_prime == "forward":
            return left_span[1] == target_index
        else:
            return right_span[0] == target_index

    def _parse_primer(self, template: str, start: int, end: int, strand: str) -> str:
        """좌표를 기반으로 프라이머 서열을 파싱하고, Reverse인 경우 역상보 변환을 수행합니다."""
        seg = template[start : end + 1].upper()
        return seg if strand == "forward" else self._rc(seg)
    
    def _build_alignment_visual(
        self,
        templates: Dict[str, str], 
        left_span: Tuple[int, int], 
        right_span: Tuple[int, int], 
        fixed_prime: str
    ) -> List[str]:
        """[AS-PCR 전용] 4가지 템플릿과 Mismatch가 반영된 모든 프라이머 세트를 시각화합니다."""
        lines = []
        
        wt_full = templates.get("wt", "")
        alt_full = templates.get("alt", "")
        wt_mm_full = templates.get("wt_mm", "")
        alt_mm_full = templates.get("alt_mm", "")

        if not wt_full or not alt_full:
            return ["ERROR: Missing Reference or Target(ALT) template sequence."]

        l_start, l_end = left_span
        r_start, r_end = right_span

        # 1. 앰플리콘 서열 추출
        wt_amp = wt_full[l_start:r_end + 1] if wt_full else ""
        alt_amp = alt_full[l_start:r_end + 1] if alt_full else ""
        wt_mm_amp = wt_mm_full[l_start:r_end + 1] if wt_mm_full else ""
        alt_mm_amp = alt_mm_full[l_start:r_end + 1] if alt_mm_full else ""

        lines.append(f"REF_SEQ            : {wt_amp}")
        lines.append(f"REF_SEQ(MISMATCH)  : {wt_mm_amp} (Pos: -{self.mismatch_pos}, Int: {self.mismatch_intensity})")
        lines.append(f"AMPLICON           : {alt_amp} (Target ALT)")
        lines.append(f"AMPLICON(MISMATCH) : {alt_mm_amp} (Target ALT)")
        lines.append("-" * (21 + len(wt_amp))) 
        
        # 2. 프라이머 파싱 헬퍼 함수
        def get_fwd(tpl): return tpl[l_start:l_end + 1] if tpl else ""
        def get_rev(tpl): return self._rc(tpl[r_start:r_end + 1]) if tpl else ""

        fwd_wt, fwd_alt = get_fwd(wt_full), get_fwd(alt_full)
        fwd_wt_mm, fwd_alt_mm = get_fwd(wt_mm_full), get_fwd(alt_mm_full)
        
        rev_wt, rev_alt = get_rev(wt_full), get_rev(alt_full)
        rev_wt_mm, rev_alt_mm = get_rev(wt_mm_full), get_rev(alt_mm_full)

        pad_rev = " " * (len(wt_amp) - len(rev_wt)) if wt_amp and rev_wt else ""

        # 3. 고정 방향(fixed_prime)에 따라 변이 프라이머 4종 출력
        if fixed_prime == "forward":
            lines.append(f"FORWARD(WT)        : {fwd_wt}")
            lines.append(f"FORWARD(WT_MM)     : {fwd_wt_mm}")
            lines.append(f"FORWARD(ALT)       : {fwd_alt}")
            lines.append(f"FORWARD(ALT_MM)    : {fwd_alt_mm}")
            lines.append(f"REVERSE(COMMON)    : {pad_rev}{self._rc(rev_wt)}")
        else:
            lines.append(f"FORWARD(COMMON)    : {fwd_wt}")
            lines.append(f"REVERSE(WT)        : {pad_rev}{self._rc(rev_wt)}")
            lines.append(f"REVERSE(WT_MM)     : {pad_rev}{self._rc(rev_wt_mm)}")
            lines.append(f"REVERSE(ALT)       : {pad_rev}{self._rc(rev_alt)}")
            lines.append(f"REVERSE(ALT_MM)    : {pad_rev}{self._rc(rev_alt_mm)}")
        
        return lines

    def _configure_force_anchor(self) -> None:
        """[AS-PCR 전용] 변이 위치(target_index)에 프라이머의 3' 말단을 강제 고정합니다."""
        self.seq_args.pop("SEQUENCE_TARGET", None)

        if self.fixed_prime == "forward":
            self.seq_args.update({"SEQUENCE_FORCE_LEFT_END": self.target_index})
            self.seq_args.pop("SEQUENCE_FORCE_LEFT_START", None)
            self.seq_args.pop("SEQUENCE_FORCE_RIGHT_START", None)
            self.seq_args.pop("SEQUENCE_FORCE_RIGHT_END", None)
        else:
            self.seq_args.update({"SEQUENCE_FORCE_RIGHT_END": self.target_index})
            self.seq_args.pop("SEQUENCE_FORCE_RIGHT_START", None)
            self.seq_args.pop("SEQUENCE_FORCE_LEFT_END", None)
            self.seq_args.pop("SEQUENCE_FORCE_LEFT_START", None)

        primer_kw = self.config.pcr_params.primer_kwargs
        n_primers = getattr(primer_kw, "n_candidates", 5)
        self.global_args.update({
            "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,"PRIMER_PICK_INTERNAL_PRIMER": 0,
            "PRIMER_NUM_RETURN": int(n_primers),
            "PRIMER_PRODUCT_SIZE_RANGE": [[getattr(primer_kw, 'min_amplicon_length', 60), getattr(primer_kw, 'max_amplicon_length', 150)]],
            "PRIMER_EXPLAIN_FLAG": 1,
        })

    # =========================================================================
    # Main Design Logic
    # =========================================================================
    def design(self) -> BaseDesignOutput:
        self._configure_force_anchor()

        res = primer3.bindings.design_primers(self.seq_args, self.global_args)
        num = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        
        flat_amplicons: List[Amplicon] = []

        if num == 0:
            error_msg = f"Primer3 Failed: L({res.get('PRIMER_LEFT_EXPLAIN', '')}), R({res.get('PRIMER_RIGHT_EXPLAIN', '')})"
            return BaseDesignOutput(status="fail", error_msg=error_msg, amplicons=[])

        chrom = self.input.reference_name
        offset = getattr(self.input, "template_genomic_start", 0)
        ref_tpl = getattr(self.input, "reference_sequence", "")
        opt_tm = float(self.global_args.get("PRIMER_OPT_TM", 58.0))

        for i in range(num):
            left_pos = res.get(f"PRIMER_LEFT_{i}")
            right_pos = res.get(f"PRIMER_RIGHT_{i}")
            if not left_pos or not right_pos: continue

            # 1. Span 계산 로직 적용
            left_span = self._left_pos_to_span(left_pos)
            right_span = self._right_pos_to_span(right_pos)
            
            # 2. 앵커 무결성 재확인
            if not self._validate_anchor(self.fixed_prime, left_span, right_span, self.target_index):
                continue
                
            # 3. WT 템플릿 파싱 무결성 검증
            wt_tpl = self.templates.get("wt", "")
            if wt_tpl and ref_tpl:
                if wt_tpl[left_span[0]:right_span[1]+1] != ref_tpl[left_span[0]:right_span[1]+1]:
                    continue
            set_id = f"Set_{i}"
            base_pair_penalty = float(res.get(f'PRIMER_PAIR_{i}_PENALTY', 0.0))

            # 4. 템플릿별 컴포넌트 파생
            for ttype in ["wt", "alt", "wt_mm", "alt_mm"]:
                tpl_seq = self.templates.get(ttype)
                if not tpl_seq: continue
                
                # 서열 파싱 (Mismatch가 포함된 서열이 잘려 나옵니다)
                fwd_seq = self._parse_primer(tpl_seq, left_span[0], left_span[1], "forward")
                rev_seq = self._parse_primer(tpl_seq, right_span[0], right_span[1], "reverse")

                # 🔥 Mismatch 여부에 따른 열역학(Tm, GC, Dimer) 동적 재계산
                fwd = Primer.make_primer(sequence=fwd_seq, role="forward", start_index=left_span[0], end_index=left_span[1] + 1)
                rev = Primer.make_primer(sequence=rev_seq, role="reverse", start_index=right_span[0], end_index=right_span[1] + 1)
                
                # 재계산된 Tm을 바탕으로 새로운 페널티 부여 (Mismatch가 발생하면 페널티가 증가하여 결과표에 반영됨)
                fwd.penalty = abs(fwd.tm - opt_tm)
                rev.penalty = abs(rev.tm - opt_tm)
                current_pair_penalty = fwd.penalty + rev.penalty + base_pair_penalty
                
                # ==========================================================
                # 공통 메타데이터 주입
                # ==========================================================
                fwd.template_sequence = tpl_seq
                fwd.reference_template_sequence = ref_tpl
                fwd.target_start_index = self.input.target_start
                fwd.target_end_index = self.input.target_end
                
                rev.template_sequence = tpl_seq
                rev.reference_template_sequence = ref_tpl
                rev.target_start_index = self.input.target_start
                rev.target_end_index = self.input.target_end

                is_fwd_as = (self.fixed_prime == "forward")
                mm_offset = -self.mismatch_pos if self.mismatch_pos else None 

                fwd.is_allele_specific = is_fwd_as
                fwd.terminal_base = fwd_seq[-1] if is_fwd_as else None
                fwd.mismatch_base = fwd_seq[mm_offset] if (is_fwd_as and mm_offset is not None) else None
                
                rev.is_allele_specific = not is_fwd_as
                rev.terminal_base = rev_seq[-1] if not is_fwd_as else None
                # 역방향 프라이머(rev_seq)도 5'->3' 기준이므로 음수 인덱스로 3' 말단을 정확히 타겟팅합니다.
                rev.mismatch_base = rev_seq[mm_offset] if (not is_fwd_as and mm_offset is not None) else None
            
                fwd.region = GenomicRegion(chrom=chrom, start=offset + left_span[0] + 1, end=offset + left_span[1] + 1, strand="+")
                rev.region = GenomicRegion(chrom=chrom, start=offset + right_span[0] + 1, end=offset + right_span[1] + 1, strand="-")
                
                alignment_view = self._build_alignment_visual(
                    self.templates, left_span, right_span, self.fixed_prime
                )
                amp = Amplicon(
                    id=f"{self.input.name}_{set_id}_{ttype}",
                    forward=fwd, reverse=rev, probe=None, 
                    template_sequence=tpl_seq,
                    target_start_index=self.input.target_start, 
                    target_end_index=self.input.target_end,
                    pair_penalty=current_pair_penalty,
                    allele_type=ttype, set_id=set_id,
                    fixed_prime=self.fixed_prime, alignment_visual=alignment_view
                )
                
                amp.is_qc_pass = True 
                flat_amplicons.append(amp)

        return BaseDesignOutput(status="success", amplicons=flat_amplicons, metadata={"assay": self.ASSAY_TYPE})