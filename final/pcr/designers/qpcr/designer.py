import primer3
from typing import Dict, Any, List

from ..base.designer import BasePrimerDesigner
from ..base.schema import BaseDesignInput, BaseDesignOutput
from pcr.components.primer import Primer, Probe
from pcr.components.amplicon import Amplicon
from pcr.components.region import GenomicRegion  # 🔥 Region 임포트 추가

class QPCRPrimerDesigner(BasePrimerDesigner):
    ASSAY_TYPE = "qPCR"

    def __init__(self, input_data: BaseDesignInput):
        super().__init__(input_data)
        self.probe_gap = self.config.pcr_params.probe_kwargs.probe_gap
        
        self.qc_stats = {"total": 0, "coverage": 0, "5_prime_g": 0, "poly_g": 0, "pass": 0} 

    # ------------------------------------------------------------------
    # 유틸리티: Alignment 텍스트 다이어그램 생성
    # ------------------------------------------------------------------
    @staticmethod
    def _reverse_complement(seq: str) -> str:
        return seq.translate(str.maketrans('ATGCatgcNn', 'TACGtacgNn'))[::-1]

    @staticmethod
    def _build_alignment_visual(ref_full: str, alt_full: str, fwd_seq: str, rev_seq: str, prb_seq: str) -> List[str]:
        lines = []
        if not ref_full or not alt_full:
            return ["ERROR: Missing Reference or Target(ALT) template sequence."]

        fwd_idx = alt_full.find(fwd_seq)
        rev_rc = QPCRPrimerDesigner._reverse_complement(rev_seq)
        rev_idx = alt_full.find(rev_rc)
        
        if fwd_idx == -1 or rev_idx == -1:
            return ["ERROR: Primers not found in the target sequence."]
            
        amp_end = rev_idx + len(rev_rc)
        
        ref_amp = ref_full[fwd_idx:amp_end]
        alt_amp = alt_full[fwd_idx:amp_end]
        
        prb_idx = alt_full.find(prb_seq)
        prb_strand = "+"
        if prb_idx == -1:
            prb_idx = alt_full.find(QPCRPrimerDesigner._reverse_complement(prb_seq))
            prb_strand = "-"

        lines.append(f"REF_SEQ  : {ref_amp}")
        lines.append(f"AMPLICON : {alt_amp} (Target ALT)")
        lines.append(f"FORWARD  : {fwd_seq}")
        
        if prb_idx != -1 and fwd_idx <= prb_idx < amp_end:
            pad = " " * (prb_idx - fwd_idx)
            lines.append(f"PROBE({prb_strand}) : {pad}{prb_seq}")
        else:
            lines.append(f"PROBE    : [NOT ALIGNED]")
            
        pad_rev = " " * (len(alt_amp) - len(rev_rc))
        lines.append(f"REVERSE  : {pad_rev}{rev_rc} (RC of {rev_seq})")
        
        return lines

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
            
            # 🔥 [속도 최적화] Primer 설계 전 중복 Probe 사전 제거 (BLAST 병목 방지)
            unique_probes = []
            seen_probe_seqs = set()
            for p in raw_probes:
                if p.sequence not in seen_probe_seqs:
                    unique_probes.append(p)
                    seen_probe_seqs.add(p.sequence)

            final_amplicons = []
            processed = 0

            # 중복이 제거된 unique_probes만 순회합니다.
            for probe_idx, probe_obj in enumerate(unique_probes):
                if processed >= limit:
                    break
                if not self._is_valid_probe(probe_obj):
                    continue
                
                amps = self._design_primers_for_probe(probe_obj, probe_idx)

                if amps:
                    final_amplicons.extend(amps)
                    processed += 1

            print(self.qc_stats)
            if not final_amplicons:
                return BaseDesignOutput(status="no_valid_pairs")
            return BaseDesignOutput(
                amplicons=final_amplicons,
                status="success",
                log_messages=[f"Raw Probes: {len(raw_probes)}, Unique: {len(unique_probes)}", f"Final Amplicons: {len(final_amplicons)}"]
            )

        except Exception as e:
            import traceback; traceback.print_exc()
            return BaseDesignOutput(status="error", error_msg=str(e))

    def _design_probes_only(self) -> List[Probe]:
        primer_opt_tm = self.config.pcr_params.primer_kwargs.opt_tm
        probe_global_args = self.config.pcr_params.probe_kwargs.to_global_args(primer_opt_tm=primer_opt_tm)
        
        run_args = self.global_args.copy()
        run_args.update(probe_global_args)
        
        n_candidates = self.config.pcr_params.probe_kwargs.n_candidates

        run_args.update({
            "PRIMER_PICK_LEFT_PRIMER":     0,
            "PRIMER_PICK_RIGHT_PRIMER":    0,
            "PRIMER_PICK_INTERNAL_OLIGO":  1,
            "PRIMER_INTERNAL_NUM_RETURN":  n_candidates,
        })

        target_len = len(self.input.alt) if getattr(self.input, "alt", None) else 1
        probe_seq_args = self.config.pcr_params.probe_kwargs.to_seq_args(
            template_seq=self.input.template_sequence,
            target_start=self.input.target_start,
            target_len=target_len
        )
        
        current_seq_args = self.seq_args.copy()
        current_seq_args.update(probe_seq_args)

        res  = primer3.bindings.design_primers(current_seq_args, run_args)
        num  = res.get("PRIMER_INTERNAL_NUM_RETURNED", 0)
        return [p for i in range(num) if (p := Probe.from_primer3(res, i))]

    def _is_valid_probe(self, probe: Probe) -> bool:
        self.qc_stats["total"] += 1
        seq = probe.sequence.upper()
        criteria = getattr(self.config.qc_criteria, "probe", None)
        if not criteria: 
            self.qc_stats["pass"] += 1
            return True
            
        p_start, p_end = probe.start_index, probe.start_index + len(seq)
        
        if not (p_start <= self.input.target_start and p_end >= self.input.target_end):
            self.qc_stats["coverage"] += 1
            return False
            
        if criteria.avoid_5_prime_g and seq.startswith("G"):
            self.qc_stats["5_prime_g"] += 1
            return False
            
        if "G" * (criteria.max_probe_poly_g) in seq:
            self.qc_stats["poly_g"] += 1
            return False
            
        self.qc_stats["pass"] += 1
        return True
    
    def _design_primers_for_probe(self, probe: Probe, probe_idx: int) -> List[Amplicon]:
        primer_kw = self.config.pcr_params.primer_kwargs
        probe_kw = self.config.pcr_params.probe_kwargs
        
        min_diff = probe_kw.min_tm_diff
        max_diff = probe_kw.max_tm_diff
        
        primer_max_tm = probe.tm - min_diff
        primer_min_tm = probe.tm - max_diff
        primer_opt_tm = (primer_max_tm + primer_min_tm) / 2.0

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
            
        res = primer3.bindings.design_primers(seq_args, run_args)
        num = res.get("PRIMER_PAIR_NUM_RETURNED", 0)

        amplicons = []
        for i in range(num):
            fwd, rev, prb = Primer.from_primer3(res, i, "LEFT"), Primer.from_primer3(res, i, "RIGHT"), Probe.from_primer3(res, i)
            if fwd and rev and prb:
                
                # 🔥 [핵심 반영] Factory에서 넘겨받은 절대 좌표(Offset) 가져오기
                chrom = self.input.reference_name
                offset = getattr(self.input, "template_genomic_start", 0)
                
                # 🔥 [핵심 반영] 각 Oligo에 방향성을 고려한 절대 게놈 좌표(GenomicRegion) 할당
                if offset > 0:
                    fwd.region = GenomicRegion(
                        chrom=chrom, 
                        start=offset + fwd.start_index + 1, 
                        end=offset + getattr(fwd, 'end_index', fwd.start_index + len(fwd.sequence)), 
                        strand="+"
                    )
                    rev.region = GenomicRegion(
                        chrom=chrom, 
                        start=offset + rev.start_index + 1, 
                        end=offset + getattr(rev, 'end_index', rev.start_index + len(rev.sequence)), 
                        strand="-"
                    )
                    prb.region = GenomicRegion(
                        chrom=chrom, 
                        start=offset + prb.start_index + 1, 
                        end=offset + getattr(prb, 'end_index', prb.start_index + len(prb.sequence)), 
                        strand="+"
                    )

                alt_seq = self.input.template_sequence
                ref_seq = getattr(self.input, "reference_sequence", "") or ""
                
                # Alignment 다이어그램 생성
                alignment_view = []
                if ref_seq and alt_seq:
                    alignment_view = self._build_alignment_visual(
                        ref_full=ref_seq, alt_full=alt_seq, 
                        fwd_seq=fwd.sequence, rev_seq=rev.sequence, prb_seq=prb.sequence
                    )

                amp = Amplicon(
                    id=f"{self.input.name}_P{probe_idx}_{i}",
                    forward=fwd, reverse=rev, probe=prb,
                    template_sequence=alt_seq,
                    reference_sequence=ref_seq,
                    target_start_index=self.input.target_start, 
                    target_end_index=self.input.target_end,
                    reference_id=chrom,
                    pair_penalty=float(res.get(f'PRIMER_PAIR_{i}_PENALTY', 0.0)) + probe.penalty,
                )
                amp.is_qc_pass = True
                setattr(amp, "alignment_visual", alignment_view) 
                
                amplicons.append(amp)

        return amplicons