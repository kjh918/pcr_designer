"""
pcr/qc/blast.py
NCBI BLAST+ 기반 특이성(Off-target) 검증 도구
(자체 In-silico PCR 좌표 계산 로직 및 Unified Amplicon Alignment 포함)
"""
import subprocess
import tempfile
import os
from typing import List, Dict, Any

try:
    import pysam
except ImportError:
    pysam = None

# 🔥 AmpliconQCStatus 임포트 추가
from pcr.config.schema.qc import BlastHit, OffTargetAmplicon, AmpliconQCStatus
from pcr.config.schema.app import PipelineConfig
from pcr.components.amplicon import Amplicon

class BlastSpecificityChecker:
    
    OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.criteria = config.qc_criteria
        self.blastn_path = config.system.blastn_path
        
        ref_key = getattr(config.pcr_params, "reference_name", "hg38")
        if ref_key not in config.references:
            ref_key = list(config.references.keys())[0]
            
        ref_config = config.references[ref_key]
        self.blast_db_path = ref_config.blast_db_path
        self.ref_path = ref_config.fasta_path
        
        self.ref_fasta = None
        self._load_ref_genome()

    def _load_ref_genome(self):
        if not pysam: 
            return
        if self.ref_path and os.path.exists(self.ref_path):
            try:
                self.ref_fasta = pysam.FastaFile(self.ref_path)
            except Exception as e:
                print(f"⚠️ Warning: Failed to load Reference FASTA: {e}")

    def __del__(self):
        if hasattr(self, 'ref_fasta') and self.ref_fasta:
            try: self.ref_fasta.close()
            except: pass

    def run(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        if not amplicons: return []
        
        fasta_content = []
        safe_id_map = {} 
        
        for i, amp in enumerate(amplicons):
            safe_id = f"AMP_{i}"
            safe_id_map[safe_id] = amp.id
            
            fasta_content.append(f">{safe_id}|Fwd\n{amp.forward.sequence}")
            fasta_content.append(f">{safe_id}|Rev\n{amp.reverse.sequence}")
            if amp.probe:
                fasta_content.append(f">{safe_id}|Probe\n{amp.probe.sequence}")

        try:
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
                tmp.write("\n".join(fasta_content))
                tmp_path = tmp.name

            cmd = [
                self.blastn_path, "-task", "blastn-short",
                "-db", self.blast_db_path, "-outfmt", self.OUTFMT,
                "-num_alignments", str(self.criteria.blast_max_alignments),
                "-num_threads", "4", "-query", tmp_path
            ]
            res = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if os.path.exists(tmp_path): os.remove(tmp_path)
            
            if not res.stdout.strip():
                print("⚠️ [BLAST Warning] No hits found! (서열이 유전체에 존재하지 않음)")

            grouped_hits = self._parse_and_group_hits(res.stdout, safe_id_map)

            for amp in amplicons:
                amp_hits = grouped_hits.get(amp.id, [])
                
                f_hits = [h for h in amp_hits if h.qseqid == "Fwd"]
                r_hits = [h for h in amp_hits if h.qseqid == "Rev"]
                p_hits = [h for h in amp_hits if h.qseqid == "Probe"]

                result_data = self._evaluate_amplicon(amp, f_hits, r_hits, p_hits)
                amp.blast_stats = result_data
                
                # 🔥 [핵심 수정] 낡은 변수 직접 덮어쓰기 방식을 버리고, 표준 qc_status 통일 규격을 사용합니다.
                if not hasattr(amp, "qc_status") or amp.qc_status is None:
                    amp.qc_status = AmpliconQCStatus()
                
                is_pass = result_data["passed"]
                messages = [] if is_pass else [result_data.get("reason", "Failed: Specificity issue.")]
                
                amp.qc_status.add_result(
                    module_name="blast",
                    is_pass=is_pass,
                    messages=messages,
                    metrics={
                        "target_count": result_data.get("target_count", 0),
                        "off_target_count": result_data.get("off_target_count", 0),
                        "noise_count": result_data.get("noise_count", 0)
                    }
                )
                
                # 하위 호환성을 위해 qc_status 상태를 기본 속성에 동기화
                amp.is_qc_pass = amp.qc_status.is_pass
                amp.qc_log = " | ".join(amp.qc_status.fail_reasons)

        except Exception as e:
            print(f"❌ BLAST Execution Failed: {e}")
            for amp in amplicons:
                if not hasattr(amp, "qc_status") or amp.qc_status is None:
                    amp.qc_status = AmpliconQCStatus()
                
                amp.qc_status.add_result(
                    module_name="blast",
                    is_pass=False,
                    messages=[f"BLAST Error: {str(e)}"]
                )
                amp.is_qc_pass = amp.qc_status.is_pass
                amp.qc_log = " | ".join(amp.qc_status.fail_reasons)

        return amplicons

    def _parse_and_group_hits(self, stdout: str, safe_id_map: dict) -> Dict[str, List[BlastHit]]:
        grouped = {}
        for line in stdout.strip().splitlines():
            cols = line.split("\t")
            if len(cols) < 14: continue
            
            raw_qseqid = cols[0]
            if "|" not in raw_qseqid: continue
            
            safe_amp_id, oligo_type = raw_qseqid.split("|", 1)
            original_amp_id = safe_id_map.get(safe_amp_id)
            if not original_amp_id: continue

            try:
                hit = BlastHit(
                    qseqid=oligo_type, 
                    sseqid=cols[1], 
                    pident=float(cols[2]), 
                    length=int(cols[3]), 
                    qstart=int(cols[6]), 
                    qend=int(cols[7]),
                    sstart=int(cols[8]), 
                    send=int(cols[9]), 
                    qseq=cols[12], 
                    sseq=cols[13]
                )

                min_id = getattr(self.criteria, "min_identity", getattr(self.criteria, "min_identity_threshold", 90.0))
                if oligo_type == "Probe":
                    min_id = getattr(self.criteria, "probe_min_identity", min_id)

                if hit.pident >= min_id:
                    if original_amp_id not in grouped: grouped[original_amp_id] = []
                    grouped[original_amp_id].append(hit)
                    
            except Exception as e:
                print(f"⚠️ [DEBUG] BlastHit Parsing Error: {e} | Line: {line}")
                continue
                
        return grouped

    def _find_amplicons(self, f_hits: List[BlastHit], r_hits: List[BlastHit]) -> List[OffTargetAmplicon]:
        amplicons = []
        for fh in f_hits:
            for rh in r_hits:
                if fh.sseqid != rh.sseqid: continue
                
                start = min(fh.genomic_start, rh.genomic_start)
                end = max(fh.genomic_end, rh.genomic_end)
                size = end - start
                
                if size <= max(500, getattr(self.criteria, "max_amp_size", 300)):
                    amplicons.append(OffTargetAmplicon(
                        chrom=fh.sseqid, start=start, end=end, 
                        product_size=size, fwd_hit=fh, rev_hit=rh,
                        is_target=False, probe_binds=False
                    ))
        return amplicons

    def _evaluate_amplicon(self, amplicon: Amplicon, f_hits: List[BlastHit], r_hits: List[BlastHit], p_hits: List[BlastHit]) -> Dict[str, Any]:
        candidates = self._find_amplicons(f_hits, r_hits)
        candidates.sort(key=lambda c: c.fwd_hit.pident + c.rev_hit.pident, reverse=True)
        
        target_signals, off_target_signals, amplification_only = [], [], []
        f_len = len(amplicon.forward.sequence)
        r_len = len(amplicon.reverse.sequence)
        has_probe = amplicon.probe is not None

        min_cov = getattr(self.criteria, "min_query_coverage", 0.8)
        min_id = getattr(self.criteria, "min_identity_threshold", 90.0)
        tolerance = getattr(self.criteria, "end_match_tolerance", 1)

        target_found = False

        for cand in candidates:
            is_convergent = False
            is_same_strand = (cand.fwd_hit.strand == cand.rev_hit.strand)
            
            if cand.fwd_hit.strand == "+" and cand.rev_hit.strand == "-": 
                if cand.fwd_hit.genomic_start <= cand.rev_hit.genomic_end:
                    is_convergent = True
            elif cand.fwd_hit.strand == "-" and cand.rev_hit.strand == "+": 
                if cand.rev_hit.genomic_start <= cand.fwd_hit.genomic_end:
                    is_convergent = True

            f_cov = cand.fwd_hit.length / f_len
            f_3end = (cand.fwd_hit.qend >= f_len - tolerance)
            f_valid = (f_cov >= min_cov) and f_3end and (cand.fwd_hit.pident >= min_id)

            r_cov = cand.rev_hit.length / r_len
            r_3end = (cand.rev_hit.qend >= r_len - tolerance)
            r_valid = (r_cov >= min_cov) and r_3end and (cand.rev_hit.pident >= min_id)

            if not (f_valid and r_valid):
                continue
                
            cand.probe_binds = self._check_probe_binding(cand, p_hits) if has_probe else True
            
            if not target_found and cand.probe_binds:
                cand.is_target = True
                target_found = True
            
            binding_report = self._create_binding_report(cand, amplicon, p_hits if has_probe and cand.probe_binds else [])
            
            if not (self.criteria.min_amp_size <= cand.product_size <= self.criteria.max_amp_size):
                binding_report["size_error"] = True
            if not is_convergent:
                binding_report["orientation_error"] = True
                binding_report["is_same_strand"] = is_same_strand

            if cand.is_target:
                target_signals.append(binding_report)
            else:
                if is_convergent:
                    if cand.probe_binds:
                        off_target_signals.append(binding_report)
                    else:
                        amplification_only.append(binding_report)

        is_passed = True
        reason = ""
        
        if len(target_signals) == 0:
            is_passed = False
            reason = "Failed: Main target amplification failed (Check Primer/Probe binding)."
        elif len(off_target_signals) > 0 or len(amplification_only) > 0:
            is_passed = False
            reason = f"Failed: Off-target detected ({len(off_target_signals)} probe-binds, {len(amplification_only)} amp-only)."
        else:
            target_report = target_signals[0]
            if target_report.get("orientation_error"):
                is_passed = False
                if target_report.get("is_same_strand"):
                    reason = "Failed: Both primers bind to the same strand. Did you forget to Reverse Complement the Reverse Primer?"
                else:
                    reason = "Failed: Primers are facing outward (Divergent). PCR will not amplify."
            elif target_report.get("size_error"):
                is_passed = False
                reason = f"Failed: Amplicon size ({target_report['product_size']}bp) is out of range ({self.criteria.min_amp_size}-{self.criteria.max_amp_size}bp)."

        return {
            "passed": is_passed, 
            "reason": reason, 
            "target_count": len(target_signals),
            "off_target_count": len(off_target_signals), 
            "noise_count": len(amplification_only),
            "target_signals": target_signals, 
            "off_target_signals": off_target_signals,
            "amplification_only": amplification_only, 
            "total_signal_count": len(target_signals) + len(off_target_signals)
        }

    def _rc(self, seq: str) -> str:
        """역상보 서열 생성 유틸"""
        return seq.translate(str.maketrans('ATGCatgcNn', 'TACGtacgNn'))[::-1]

    def _fetch_sequence_from_genome(self, chrom: str, start: int, end: int) -> str:
        """FASTA로부터 서열을 가져옵니다."""
        if not self.ref_fasta: return ""
        try:
            return self.ref_fasta.fetch(chrom, start, end).upper()
        except:
            alt = chrom.replace("chr", "") if "chr" in chrom else f"chr{chrom}"
            try: return self.ref_fasta.fetch(alt, start, end).upper()
            except: return ""

    def _create_binding_report(self, cand: OffTargetAmplicon, amplicon: Amplicon, p_hits: List[BlastHit]) -> Dict[str, Any]:
        """
        프론트엔드 UI 표출(Amplicon block)에 필요한 텍스트(unified_text_block) 1개만 남기고 
        무거운 서열들은 제거된 경량화 버전
        """
        fwd_seq = amplicon.forward.sequence
        rev_seq = amplicon.reverse.sequence
        prb_seq = amplicon.probe.sequence if amplicon.probe else ""

        def get_bounds(hit, seq):
            if hit.strand == "+":
                return hit.genomic_start - (hit.qstart - 1), hit.genomic_end + (len(seq) - hit.qend)
            else:
                return hit.genomic_start - (len(seq) - hit.qend), hit.genomic_end + (hit.qstart - 1)

        f_gstart, f_gend = get_bounds(cand.fwd_hit, fwd_seq)
        r_gstart, r_gend = get_bounds(cand.rev_hit, rev_seq)

        amp_gstart = min(f_gstart, r_gstart)
        amp_gend = max(f_gend, r_gend)

        target_p_hit = next((ph for ph in p_hits if ph.genomic_start >= cand.start and ph.genomic_end <= cand.end), None)
        if target_p_hit and prb_seq:
            p_gstart, p_gend = get_bounds(target_p_hit, prb_seq)
            amp_gstart = min(amp_gstart, p_gstart)
            amp_gend = max(amp_gend, p_gend)

        full_ref_seq = self._fetch_sequence_from_genome(cand.chrom, amp_gstart - 1, amp_gend)
        if not full_ref_seq:
            full_ref_seq = "N" * (amp_gend - amp_gstart + 1)
        elif len(full_ref_seq) < (amp_gend - amp_gstart + 1):
            full_ref_seq = full_ref_seq.ljust((amp_gend - amp_gstart + 1), 'N')

        overall_match = [" "] * len(full_ref_seq)

        def format_oligo(hit, seq):
            h_gstart, _ = get_bounds(hit, seq)
            offset = h_gstart - amp_gstart
            
            aligned_seq = seq.upper() if hit.strand == "+" else self._rc(seq.upper())
            
            if hit.strand == "+":
                pad_left = hit.qstart - 1
                pad_right = len(seq) - hit.qend
            else:
                pad_left = len(seq) - hit.qend
                pad_right = hit.qstart - 1

            formatted_input = ""
            for i in range(len(aligned_seq)):
                is_overhang = (i < pad_left) or (i >= len(seq) - pad_right)
                ref_idx = offset + i
                
                if 0 <= ref_idx < len(full_ref_seq):
                    ref_base = full_ref_seq[ref_idx]
                    if aligned_seq[i].upper() == ref_base.upper() and ref_base.upper() != 'N':
                        formatted_input += aligned_seq[i].upper()
                        if not is_overhang:
                            overall_match[ref_idx] = "|"
                        else:
                            formatted_input = formatted_input[:-1] + aligned_seq[i].lower()
                            overall_match[ref_idx] = "."
                    else:
                        formatted_input += aligned_seq[i].lower()
                        overall_match[ref_idx] = "."
                else:
                    formatted_input += aligned_seq[i].lower()

            input_line = (" " * offset) + formatted_input
            blast_line = (" " * (offset + pad_left)) + hit.qseq
            
            return input_line, blast_line

        f_input, f_blast = format_oligo(cand.fwd_hit, fwd_seq)
        r_input, r_blast = format_oligo(cand.rev_hit, rev_seq)
        
        match_line = "".join(overall_match)

        title_prefix = "🎯 TARGET" if cand.is_target else "⚠️ NON-SPECIFIC"
        
        lines = [
            f"[{title_prefix} Amplicon Alignment | {cand.chrom}:{amp_gstart}-{amp_gend} | Size: {amp_gend - amp_gstart + 1}bp]",
            f"REF_SEQ : {full_ref_seq}",
            f"MATCH   : {match_line}",
            f"FORWARD : {f_input} ({cand.fwd_hit.strand})",
            f"F_BLAST : {f_blast}",
            f"REVERSE : {r_input} ({cand.rev_hit.strand})",
            f"R_BLAST : {r_blast}"
        ]

        if target_p_hit and prb_seq:
            p_input, p_blast = format_oligo(target_p_hit, prb_seq)
            lines.extend([
                f"PROBE   : {p_input} ({target_p_hit.strand})",
                f"P_BLAST : {p_blast}"
            ])

        f_report = {
            "label": "Forward Primer", "identity": f"{cand.fwd_hit.pident}%",
            "chrom": cand.fwd_hit.sseqid, "strand": cand.fwd_hit.strand,
            "coordinate": f"{cand.fwd_hit.sseqid}:{cand.fwd_hit.genomic_start}-{cand.fwd_hit.genomic_end} ({cand.fwd_hit.strand})"
        }
        r_report = {
            "label": "Reverse Primer", "identity": f"{cand.rev_hit.pident}%",
            "chrom": cand.rev_hit.sseqid, "strand": cand.rev_hit.strand,
            "coordinate": f"{cand.rev_hit.sseqid}:{cand.rev_hit.genomic_start}-{cand.rev_hit.genomic_end} ({cand.rev_hit.strand})"
        }
        p_report = None
        if target_p_hit:
            p_report = {
                "label": "TaqMan Probe", "identity": f"{target_p_hit.pident}%",
                "chrom": target_p_hit.sseqid, "strand": target_p_hit.strand,
                "coordinate": f"{target_p_hit.sseqid}:{target_p_hit.genomic_start}-{target_p_hit.genomic_end} ({target_p_hit.strand})"
            }

        return {
            "location": f"{cand.chrom}:{amp_gstart}-{amp_gend}",
            "product_size": amp_gend - amp_gstart + 1, 
            "is_target": cand.is_target, 
            "fwd": f_report, 
            "rev": r_report, 
            "probe": p_report,
            "unified_text_block": "\n".join(lines) 
        }

    def _check_probe_binding(self, cand: OffTargetAmplicon, p_hits: List[BlastHit]) -> bool:
        for ph in p_hits:
            if ph.sseqid != cand.chrom: continue
            if ph.genomic_start >= cand.start and ph.genomic_end <= cand.end: return True
        return False