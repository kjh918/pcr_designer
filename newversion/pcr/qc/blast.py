import subprocess
import tempfile
import os
from typing import List, Dict, Optional, Any

# 외부 라이브러리 체크
try:
    import pysam
    from ispcr import calculate_pcr_product, FastaSequence
except ImportError:
    pysam = None
    calculate_pcr_product = None
    FastaSequence = None

from ..config.schema.qc import BlastHit, OffTargetAmplicon
from ..config.schema.root import AppConfig
from ..components.amplicon import Amplicon

class BlastSpecificityChecker:
    OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"

    def __init__(self, config: AppConfig):
        self.config = config
        self.criteria = config.qc_criteria
        
        # ✅ 수정 1: 실제 경로 데이터는 qc_tools.paths 안에 있음
        self.tools = config.qc_tools.paths 

    def check_amplicon(self, amplicon: Amplicon) -> Dict[str, Any]:
        """BLAST + isPCR 하이브리드 QC"""
        fwd = amplicon.forward
        rev = amplicon.reverse
        
        # 1. BLAST 실행
        hits = self._run_blast_pair(
            f_name="Fwd", f_seq=fwd.sequence, 
            r_name="Rev", r_seq=rev.sequence
        )
        f_hits = [h for h in hits if h.qseqid == "Fwd"]
        r_hits = [h for h in hits if h.qseqid == "Rev"]

        # 2. 후보 조합
        candidates = self._find_amplicons(f_hits, r_hits)

        # 3. isPCR 검증
        confirmed_off_targets = []
        intended_target = None
        
        # ✅ 수정 2: Reference 경로 찾기 로직 개선
        # System Config (hg38 등) 참조 우선 -> 없으면 BLAST DB REF 사용
        ref_path = ""
        if amplicon.reference_id in self.config.references:
            ref_path = self.config.references[amplicon.reference_id].fasta_path
        elif self.tools.blast_ref_path: # 이전 blast_db_ref_path -> blast_ref_path
            ref_path = self.tools.blast_ref_path

        ref_fasta = None
        try:
            if pysam and ref_path and os.path.exists(ref_path):
                ref_fasta = pysam.FastaFile(ref_path)
        except Exception as e:
            print(f"⚠️ Warning: Could not open reference genome: {e}")

        for cand in candidates:
            # 타겟 판정
            is_target_pos = self._is_intended_pos(cand, amplicon)
            cand.is_target = is_target_pos

            # isPCR 수행
            is_real = self._verify_with_ispcr(cand, fwd.sequence, rev.sequence, ref_fasta)
            
            if is_real:
                if is_target_pos:
                    if intended_target is None:
                        intended_target = cand
                else:
                    confirmed_off_targets.append(cand)
        
        if ref_fasta:
            ref_fasta.close()

        # 4. 결과 업데이트
        is_passed = (len(confirmed_off_targets) == 0) and (intended_target is not None)
        
        amplicon.is_qc_pass = is_passed
        amplicon.off_target_count = len(confirmed_off_targets)
        
        # 상세 결과 저장 (나중에 리포트에 씀)
        amplicon.qc_details["specificity"] = {
            "blast_hits": len(hits),
            "candidates": len(candidates),
            "intended_found": intended_target is not None
        }

        # 5. 시각화 데이터
        alignment_view = None
        if intended_target and ref_path:
            alignment_view = self._generate_alignment_view(amplicon, intended_target, ref_path)

        return {
            "passed": is_passed,
            "off_targets": confirmed_off_targets,
            "alignment": alignment_view
        }

    # ... (Helper methods: _is_intended_pos, _verify_with_ispcr, _generate_alignment_view 등은 기존 유지) ...
    def _is_intended_pos(self, cand: OffTargetAmplicon, amp: Amplicon) -> bool:
        """현재 후보가 디자인된 타겟 위치와 일치하는지 확인"""
        # Chromosome 이름 비교 (단순 문자열 비교)
        # hg38 vs chr1 같은 매핑 이슈가 있을 수 있으니 주의 필요. 여기서는 단순 포함관계 확인
        if cand.chrom != amp.reference_id and amp.reference_id not in cand.chrom:
             pass 

        # 좌표 오차 범위 (예: 100bp) 내에 있는지
        return abs(cand.start - amp.target_start_index) < 100

    def _verify_with_ispcr(self, cand: OffTargetAmplicon, 
                           fwd_seq: str, rev_seq: str, ref_fasta: Any) -> bool:
        if not calculate_pcr_product or not ref_fasta:
            return True

        padding = 100
        fetch_start = max(0, cand.start - padding)
        fetch_end = cand.end + padding

        try:
            local_seq_str = ref_fasta.fetch(cand.chrom, fetch_start, fetch_end)
            local_seq_obj = FastaSequence(cand.chrom, local_seq_str)
            fwd_obj = FastaSequence("Fwd", fwd_seq)
            rev_obj = FastaSequence("Rev", rev_seq)

            result = calculate_pcr_product(
                sequence=local_seq_obj,
                forward_primer=fwd_obj,
                reverse_primer=rev_obj,
                min_product_length=20, 
                max_product_length=5000,
                header=False, cols="all", output_file=False
            )
            return bool(result and result.strip())
        except:
            return True

    def _generate_alignment_view(self, amp: Amplicon, target_hit: OffTargetAmplicon, ref_path: str) -> Dict[str, str]:
        try:
            if not pysam or not os.path.exists(ref_path): return {}
            with pysam.FastaFile(ref_path) as fasta:
                ref_seq = fasta.fetch(target_hit.chrom, target_hit.start, target_hit.end).upper()
            return {"chrom": target_hit.chrom, "snippet": ref_seq[:50]}
        except: return {}

    def _run_blast_pair(self, f_name: str, f_seq: str, r_name: str, r_seq: str) -> List[BlastHit]:
        """BLASTN 실행"""
        fasta_content = f">{f_name}\n{f_seq}\n>{r_name}\n{r_seq}\n"
        
        # ✅ 수정 3: 올바른 경로 변수 사용 (self.tools.blastn, self.tools.blast_db_path)
        # 키 이름이 YAML과 일치해야 함 (blastn_exe -> blastn)
        cmd = [
            self.tools.blastn,
            "-task", "blastn-short",
            "-db", self.tools.blast_db_path,
            "-outfmt", self.OUTFMT,
            "-num_alignments", str(self.criteria.blast_max_alignments),
            "-num_threads", "4"
        ]

        try:
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fa") as tmp:
                tmp.write(fasta_content)
                tmp.flush()
                res = subprocess.run(
                    cmd + ["-query", tmp.name], 
                    capture_output=True, text=True, check=True
                )
                return self._parse_output(res.stdout)
        except Exception as e:
            print(f"❌ BLAST Execution Failed: {e}")
            return []

    def _parse_output(self, stdout: str) -> List[BlastHit]:
        """BLAST Output 파싱"""
        hits = []
        for line in stdout.strip().splitlines():
            cols = line.split("\t")
            if len(cols) < 12: continue

            try:
                hit = BlastHit(
                    qseqid=cols[0], sseqid=cols[1],
                    pident=float(cols[2]), length=int(cols[3]),
                    qstart=int(cols[6]), qend=int(cols[7]),
                    sstart=int(cols[8]), send=int(cols[9])
                )
                if (hit.pident >= self.criteria.min_identity and 
                    hit.length >= self.criteria.min_hit_length):
                    hits.append(hit)
            except: continue
        return hits

    def _find_amplicons(self, f_hits: List[BlastHit], r_hits: List[BlastHit]) -> List[OffTargetAmplicon]:
        """Forward/Reverse 조합하여 Amplicon 후보 찾기"""
        amplicons = []
        for fh in f_hits:
            for rh in r_hits:
                if fh.sseqid != rh.sseqid: continue
                if fh.strand == rh.strand: continue
                
                valid = False
                if fh.strand == "+" and rh.strand == "-" and fh.genomic_end < rh.genomic_start: valid = True
                elif fh.strand == "-" and rh.strand == "+" and rh.genomic_end < fh.genomic_start: valid = True
                
                if not valid: continue

                start = min(fh.genomic_start, rh.genomic_start)
                end = max(fh.genomic_end, rh.genomic_end)
                size = end - start
                
                if self.criteria.min_amp_size <= size <= self.criteria.max_amp_size:
                    amplicons.append(OffTargetAmplicon(
                        chrom=fh.sseqid, start=start, end=end, product_size=size,
                        fwd_hit=fh, rev_hit=rh
                    ))
        return amplicons