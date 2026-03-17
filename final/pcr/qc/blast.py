"""
pcr/qc/blast.py
NCBI BLAST+ 기반 특이성(Off-target) 검증 도구
(자체 In-silico PCR 좌표 계산 로직 포함)
"""
import subprocess
import tempfile
import os
from typing import List, Dict, Any

try:
    import pysam
except ImportError:
    pysam = None

from pcr.config.schema.qc import BlastHit, OffTargetAmplicon
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
                amp.is_qc_pass = result_data["passed"]
                
                if amp.is_qc_pass:
                    amp.qc_log = "Passed: Specific target confirmed."
                else:
                    amp.qc_log = result_data.get("reason", "Failed: Specificity issue.")

        except Exception as e:
            print(f"❌ BLAST Execution Failed: {e}")
            for amp in amplicons:
                amp.is_qc_pass = False
                amp.qc_log = f"Error: BLAST analysis failed ({str(e)})"

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
        """
        🔥 [스마트 클러스터링]
        가닥(Strand) 방향을 먼저 픽스하지 않고, 거리가 가까우면 무조건 앰플리콘 후보로 묶어줍니다.
        나중에 방향성(Orientation)을 검증하여 사용자에게 친절한 에러를 제공합니다.
        """
        amplicons = []
        for fh in f_hits:
            for rh in r_hits:
                if fh.sseqid != rh.sseqid: continue
                
                # 단순히 거리가 가까운 것들을 모두 추출
                start = min(fh.genomic_start, rh.genomic_start)
                end = max(fh.genomic_end, rh.genomic_end)
                size = end - start
                
                if size <= max(5000, getattr(self.criteria, "max_amp_size", 300)):
                    amplicons.append(OffTargetAmplicon(
                        chrom=fh.sseqid, 
                        start=start, 
                        end=end, 
                        product_size=size, 
                        fwd_hit=fh, 
                        rev_hit=rh,
                        is_target=False,
                        probe_binds=False
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
            # 1. 방향성(Orientation) 판별: 마주보고 있는가(Convergent)?
            is_convergent = False
            is_same_strand = (cand.fwd_hit.strand == cand.rev_hit.strand)
            
            if cand.fwd_hit.strand == "+" and cand.rev_hit.strand == "-": 
                if cand.fwd_hit.genomic_start <= cand.rev_hit.genomic_end:
                    is_convergent = True
            elif cand.fwd_hit.strand == "-" and cand.rev_hit.strand == "+": 
                if cand.rev_hit.genomic_start <= cand.fwd_hit.genomic_end:
                    is_convergent = True

            # 2. 프라이머 결합 검증
            f_cov = cand.fwd_hit.length / f_len
            f_3end = (cand.fwd_hit.qend >= f_len - tolerance)
            f_valid = (f_cov >= min_cov) and f_3end and (cand.fwd_hit.pident >= min_id)

            r_cov = cand.rev_hit.length / r_len
            r_3end = (cand.rev_hit.qend >= r_len - tolerance)
            r_valid = (r_cov >= min_cov) and r_3end and (cand.rev_hit.pident >= min_id)

            if not (f_valid and r_valid):
                continue
                
            # 3. 프로브 결합 여부 (프로브가 없으면 True 간주)
            cand.probe_binds = self._check_probe_binding(cand, p_hits) if has_probe else True
            
            # 4. 진짜 타겟 식별
            if not target_found and cand.probe_binds:
                cand.is_target = True
                target_found = True
            
            binding_report = self._create_binding_report(cand, amplicon, p_hits if has_probe and cand.probe_binds else [])
            
            # 에러 마킹
            if not (self.criteria.min_amp_size <= cand.product_size <= self.criteria.max_amp_size):
                binding_report["size_error"] = True
            if not is_convergent:
                binding_report["orientation_error"] = True
                binding_report["is_same_strand"] = is_same_strand

            # 5. 신호 분류
            if cand.is_target:
                target_signals.append(binding_report)
            else:
                # 🔥 오프타겟의 경우 PCR 증폭이 실제로 가능한(Convergent) 형태일 때만 기록 (노이즈 방지)
                if is_convergent:
                    if cand.probe_binds:
                        off_target_signals.append(binding_report)
                    else:
                        amplification_only.append(binding_report)

        # 6. 최종 QC 판정 로직
        is_passed = True
        reason = ""
        
        if len(target_signals) == 0:
            is_passed = False
            reason = "Failed: Main target amplification failed (Check Primer/Probe binding)."
        elif len(off_target_signals) > 0 or len(amplification_only) > 0:
            is_passed = False
            reason = f"Failed: Off-target detected ({len(off_target_signals)} probe-binds, {len(amplification_only)} amp-only)."
        else:
            # 타겟은 1개인데 결함(방향 또는 사이즈)이 있는 경우 친절하게 에러 반환
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
            "passed": is_passed, "reason": reason, "target_count": len(target_signals),
            "off_target_count": len(off_target_signals), "noise_count": len(amplification_only),
            "target_signals": target_signals, "off_target_signals": off_target_signals,
            "amplification_only": amplification_only, "total_signal_count": len(target_signals) + len(off_target_signals)
        }

    def _create_binding_report(self, cand: OffTargetAmplicon, amplicon: Amplicon, p_hits: List[BlastHit]) -> Dict[str, Any]:
        # BlastHit 객체 자체를 넘겨서 게놈 좌표 및 방향성 등의 상세 정보를 추출하게 함
        f_report = self._generate_alignment_block("Forward Primer", cand.fwd_hit)
        r_report = self._generate_alignment_block("Reverse Primer", cand.rev_hit)
        
        p_report = None
        target_p_hit = next((ph for ph in p_hits if ph.genomic_start >= cand.start and ph.genomic_end <= cand.end), None)
        if target_p_hit:
            p_report = self._generate_alignment_block("TaqMan Probe", target_p_hit)

        ref_seq = self._fetch_sequence_from_genome(cand.chrom, cand.start - 1, cand.end)

        return {
            "location": f"{cand.chrom}:{cand.start}-{cand.end}",
            "product_size": cand.product_size, 
            "is_target": cand.is_target, 
            "fwd": f_report, "rev": r_report, "probe": p_report,
            "full_sequence": ref_seq
        }

    def _generate_alignment_block(self, name: str, hit: BlastHit) -> Dict[str, Any]:
        """미스매치 시각화 및 실제 게놈 좌표 상세 정보를 포함하는 딕셔너리 생성"""
        query_seq = hit.qseq
        subject_seq = hit.sseq
        match_line = "".join("|" if q == s else " " if q == "-" or s == "-" else "." for q, s in zip(query_seq, subject_seq))
        
        return {
            "label": name, 
            "query": query_seq, 
            "match": match_line,
            "subject": subject_seq, 
            "identity": f"{hit.pident}%",
            "chrom": hit.sseqid,
            "strand": hit.strand,
            "genomic_start": hit.genomic_start,
            "genomic_end": hit.genomic_end,
            "coordinate": f"{hit.sseqid}:{hit.genomic_start}-{hit.genomic_end} ({hit.strand})"
        }

    def _check_probe_binding(self, cand: OffTargetAmplicon, p_hits: List[BlastHit]) -> bool:
        for ph in p_hits:
            if ph.sseqid != cand.chrom: continue
            if ph.genomic_start >= cand.start and ph.genomic_end <= cand.end: return True
        return False
        
    def _fetch_sequence_from_genome(self, chrom: str, start: int, end: int) -> str:
        if not self.ref_fasta: return ""
        try:
            return self.ref_fasta.fetch(chrom, start, end).upper()
        except:
            alt = chrom.replace("chr", "") if "chr" in chrom else f"chr{chrom}"
            try: return self.ref_fasta.fetch(alt, start, end).upper()
            except: return ""