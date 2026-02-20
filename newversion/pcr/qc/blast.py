import subprocess
import tempfile
import os
from typing import List, Dict, Any

try:
    import pysam
except ImportError:
    pysam = None

from ..config.schema.qc import BlastHit, OffTargetAmplicon
from ..config.schema.root import AppConfig
from ..components.amplicon import Amplicon

class BlastSpecificityChecker:
    OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"

    def __init__(self, config: AppConfig):
        self.config = config
        self.criteria = config.qc_criteria
        self.tools = config.qc_tools.paths 
        
        # Reference Fasta 로드 (서열 추출용)
        self.ref_fasta = None
        self._load_ref_genome()

    def _load_ref_genome(self):
        if not pysam: return
        
        # 1. Config에 정의된 Reference 경로 확인
        ref_path = self.tools.blast_ref_path 
        if not ref_path and self.config.references:
             first_ref = next(iter(self.config.references.values()))
             ref_path = first_ref.fasta_path

        if ref_path and os.path.exists(ref_path):
            try:
                self.ref_fasta = pysam.FastaFile(ref_path)
            except Exception as e:
                print(f"⚠️ Warning: Failed to load Reference FASTA: {e}")

    def __del__(self):
        if self.ref_fasta:
            try: self.ref_fasta.close()
            except: pass

    def check_amplicon(self, amplicon: Amplicon) -> Dict[str, Any]:
        """
        [BLAST-Only Logic]
        isPCR 라이브러리 없이, BLAST 좌표 계산만으로 증폭 여부를 판단합니다.
        """
        fwd = amplicon.forward
        rev = amplicon.reverse
        probe = amplicon.probe

        # 1. BLAST 실행 (Fwd, Rev, Probe)
        hits = self._run_blast_set(fwd.sequence, rev.sequence, probe.sequence if probe else None)
        
        f_hits = [h for h in hits if h.qseqid == "Fwd"]
        r_hits = [h for h in hits if h.qseqid == "Rev"]
        p_hits = [h for h in hits if h.qseqid == "Probe"]

        # 2. 증폭 산물 후보 찾기 (여기가 곧 in-silico PCR 역할)
        # Fwd와 Rev가 마주보고 있고, 거리가 적절하면 Amplicon으로 간주
        candidates = self._find_amplicons(f_hits, r_hits)

        # 3. Probe 검증 및 분류
        signal_candidates = []      # Probe까지 붙는 후보 (형광 O)
        amplification_candidates = [] # Probe는 안 붙는 후보 (형광 X, 증폭 O)

        for cand in candidates:
            # 서열 추출 (리포팅용)
            seq_data = self._fetch_sequence_from_genome(cand.chrom, cand.start, cand.end)
            cand.reference_sequence = seq_data
            cand.template_sequence = seq_data

            # Probe 결합 확인 (좌표 기반)
            probe_binds = self._check_probe_binding(cand, p_hits)
            cand.probe_binds = probe_binds

            if probe_binds:
                signal_candidates.append(cand)
            else:
                amplification_candidates.append(cand)

        # 4. 결과 판정 (Count 기반)
        # - Signal 후보가 딱 1개여야 함 (그게 Intended Target)
        count = len(signal_candidates)
        is_passed = (count == 1)
        
        intended_target = None
        off_targets = []

        if count == 1:
            intended_target = signal_candidates[0]
            intended_target.is_target = True
        elif count > 1:
            off_targets = signal_candidates
        
        # 5. 결과 저장
        amplicon.is_qc_pass = is_passed
        amplicon.off_target_count = len(off_targets) if count > 1 else 0
        
        amplicon.qc_details["specificity"] = {
            "blast_hits_total": len(hits),
            "signal_candidates_count": count,       # 1이면 정상, >1이면 Off-target
            "amplification_only_count": len(amplification_candidates),
            "intended_found": (count == 1)
        }

        # 시각화용 데이터
        alignment_view = {}
        target_to_show = intended_target if intended_target else (off_targets[0] if off_targets else None)
        
        if target_to_show and target_to_show.reference_sequence:
             alignment_view = {
                "chrom": target_to_show.chrom,
                "start": target_to_show.start,
                "end": target_to_show.end,
                "seq_snippet": target_to_show.reference_sequence[:50]
            }

        return {
            "passed": is_passed,
            "off_targets": off_targets,
            "alignment": alignment_view
        }

    # --------------------------------------------------------------------------
    # Helper Methods
    # --------------------------------------------------------------------------
    def _fetch_sequence_from_genome(self, chrom: str, start: int, end: int) -> str:
        if not self.ref_fasta: return ""
        try:
            # pysam fetch
            return self.ref_fasta.fetch(chrom, start, end).upper()
        except KeyError:
            # chr 처리
            alt = chrom.replace("chr", "") if "chr" in chrom else f"chr{chrom}"
            try: return self.ref_fasta.fetch(alt, start, end).upper()
            except: return ""
        except: return ""

    def _run_blast_set(self, fwd_seq: str, rev_seq: str, probe_seq: str = None) -> List[BlastHit]:
        """BLASTN 실행"""
        fasta_content = f">Fwd\n{fwd_seq}\n>Rev\n{rev_seq}\n"
        if probe_seq: fasta_content += f">Probe\n{probe_seq}\n"
        
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
                tmp.write(fasta_content); tmp.flush()
                res = subprocess.run(
                    cmd + ["-query", tmp.name], 
                    capture_output=True, text=True, check=True
                )
                return self._parse_output(res.stdout)
        except Exception as e:
            print(f"❌ BLAST Execution Failed: {e}")
            return []

    def _check_probe_binding(self, cand: OffTargetAmplicon, p_hits: List[BlastHit]) -> bool:
        """Probe가 Amplicon 내부에 존재하는지 확인"""
        for ph in p_hits:
            if ph.sseqid != cand.chrom: continue
            
            # Amplicon 범위 내에 Probe가 매핑되는지 (Overlap)
            if ph.genomic_start >= cand.start and ph.genomic_end <= cand.end:
                return True
        return False

    def _find_amplicons(self, f_hits: List[BlastHit], r_hits: List[BlastHit]) -> List[OffTargetAmplicon]:
        """
        [BLAST 기반 in-silico PCR]
        Forward와 Reverse Hit의 위치 관계를 분석하여 증폭 가능한 후보를 찾습니다.
        """
        amplicons = []
        for fh in f_hits:
            for rh in r_hits:
                # 1. 같은 염색체여야 함
                if fh.sseqid != rh.sseqid: continue
                # 2. 서로 반대 스트랜드여야 함 (PCR 원리)
                if fh.strand == rh.strand: continue
                
                valid = False
                start, end = 0, 0
                
                # Case A: Fwd(+) ... Rev(-)
                if fh.strand == "+" and rh.strand == "-" and fh.genomic_end < rh.genomic_start: 
                    valid = True
                    start, end = fh.genomic_start, rh.genomic_end
                
                # Case B: Rev(+) ... Fwd(-) (Reverse Primer가 Forward 역할 하는 경우)
                elif fh.strand == "-" and rh.strand == "+" and rh.genomic_end < fh.genomic_start: 
                    valid = True
                    start, end = rh.genomic_start, fh.genomic_end
                
                if valid:
                    size = end - start
                    # 3. Product Size 체크 (설정된 범위 내)
                    if self.criteria.min_amp_size <= size <= self.criteria.max_amp_size:
                        amplicons.append(OffTargetAmplicon(
                            chrom=fh.sseqid, 
                            start=start, 
                            end=end, 
                            product_size=size, 
                            fwd_hit=fh, 
                            rev_hit=rh
                        ))
        return amplicons
    
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
                if hit.pident >= self.criteria.min_identity: 
                    hits.append(hit)
            except: continue
        return hits