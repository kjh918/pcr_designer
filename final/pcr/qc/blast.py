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

# 프로젝트 스키마 및 컴포넌트 임포트
from pcr.config.schema.qc import BlastHit, OffTargetAmplicon
from pcr.config.schema.app import PipelineConfig
from pcr.components.amplicon import Amplicon

class BlastSpecificityChecker:
    OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    def __init__(self, config: PipelineConfig):
        """
        [MODIFIED] config(PipelineConfig)에서 시스템 및 레퍼런스 정보를 동적으로 로드합니다.
        """
        self.config = config
        self.criteria = config.qc_criteria
        
        # 1. system.yaml에 정의된 blastn 실행 경로
        self.blastn_path = config.system.blastn_path
        
        # 2. references 설정에서 첫 번째 레퍼런스(기본값 hg38 등) 정보를 가져옴
        ref_key = "hg38" if "hg38" in config.references else list(config.references.keys())[0]
        ref_config = config.references[ref_key]
        
        self.blast_db_path = ref_config.blast_db_path
        self.ref_path = ref_config.fasta_path
        
        # 3. Genomic sequence 조회를 위한 pysam 핸들
        self.ref_fasta = None
        self._load_ref_genome()

    def _load_ref_genome(self):
        """Reference FASTA 파일을 로드합니다."""
        if not pysam: 
            return
        if self.ref_path and os.path.exists(self.ref_path):
            try:
                self.ref_fasta = pysam.FastaFile(self.ref_path)
            except Exception as e:
                print(f"⚠️ Warning: Failed to load Reference FASTA: {e}")

    def __del__(self):
        """객체 소멸 시 pysam 리소스 해제"""
        if hasattr(self, 'ref_fasta') and self.ref_fasta:
            try:
                self.ref_fasta.close()
            except:
                pass

    def run(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        """
        Batch 처리를 통해 모든 후보군의 특이성을 한 번에 검증합니다.
        """
        
        if not amplicons:
            return []
        blast_db = getattr(self.config.qc_criteria, "blast_db", None)
        
        if not blast_db or blast_db.lower() == "none":
            return amplicons # 검사 없이 전원 합격 처리
        
        # 1. Multi-FASTA 생성 (Query 통합)
        fasta_content = []
        for amp in amplicons:
            fasta_content.append(f">{amp.id}|Fwd\n{amp.forward.sequence}")
            fasta_content.append(f">{amp.id}|Rev\n{amp.reverse.sequence}")
            if amp.probe:
                fasta_content.append(f">{amp.id}|Probe\n{amp.probe.sequence}")

        valid_amplicons = []

        try:
            # 2. 임시 파일 생성 및 BLAST 실행
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
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

            # 3. 결과 파싱 및 그룹화
            grouped_hits = self._parse_and_group_hits(res.stdout)

            # 4. 개별 앰플리콘 별로 In-silico PCR 결과 판정
            for amp in amplicons:
                amp_hits = grouped_hits.get(amp.id, [])
                
                f_hits = [h for h in amp_hits if h.qseqid == "Fwd"]
                r_hits = [h for h in amp_hits if h.qseqid == "Rev"]
                p_hits = [h for h in amp_hits if h.qseqid == "Probe"]

                result_data = self._evaluate_amplicon(amp, f_hits, r_hits, p_hits)
                amp.blast_stats = result_data
                
                if result_data["passed"]:
                    valid_amplicons.append(amp)
                else:
                    amp.is_qc_pass = False
                    amp.qc_log = "Failed BLAST: Non-specific amplification detected."

        except Exception as e:
            print(f"❌ BLAST Execution Failed: {e}")

        return valid_amplicons

    def _evaluate_amplicon(self, amplicon: Amplicon, f_hits: List[BlastHit], r_hits: List[BlastHit], p_hits: List[BlastHit]) -> Dict[str, Any]:
        """In-silico PCR 판별 및 오프타겟 유무 검사"""
        candidates = self._find_amplicons(f_hits, r_hits)

        signal_candidates = []
        amplification_only_candidates = []

        for cand in candidates:
            # 게놈에서 서열 추출 (pysam 활용)
            cand.reference_sequence = self._fetch_sequence_from_genome(cand.chrom, cand.start, cand.end)
            
            # TaqMan Probe가 있을 경우 바인딩 여부 확인
            probe_binds = self._check_probe_binding(cand, p_hits) if amplicon.probe else False
            
            if not amplicon.probe:
                signal_candidates.append(cand) # SYBR Green 방식 등
            else:
                if probe_binds:
                    signal_candidates.append(cand)
                else:
                    amplification_only_candidates.append(cand)

        # 판정: 시그널이 정확히 1개(의도한 타겟)만 나와야 통과
        count = len(signal_candidates)
        is_passed = (count == 1)
        
        alignment_view = {}
        target = signal_candidates[0] if count >= 1 else None
        if target and target.reference_sequence:
             alignment_view = {
                "chrom": target.chrom, "start": target.start, "end": target.end,
                "seq_snippet": target.reference_sequence[:50]
            }

        return {
            "passed": is_passed,
            "signal_candidates_count": count,
            "amplification_only_count": len(amplification_only_candidates),
            "alignment": alignment_view
        }

    def _parse_and_group_hits(self, stdout: str) -> Dict[str, List[BlastHit]]:
        grouped = {}
        for line in stdout.strip().splitlines():
            cols = line.split("\t")
            if len(cols) < 12: continue
            
            raw_qseqid = cols[0]
            if "|" not in raw_qseqid: continue
            amp_id, oligo_type = raw_qseqid.split("|", 1)

            try:
                sstart, send = int(cols[8]), int(cols[9])
                hit = BlastHit(
                    qseqid=oligo_type, sseqid=cols[1], pident=float(cols[2]), 
                    length=int(cols[3]), qstart=int(cols[6]), qend=int(cols[7]),
                    sstart=sstart, send=send
                )
                
                if hit.pident >= self.criteria.min_identity:
                    if amp_id not in grouped: grouped[amp_id] = []
                    grouped[amp_id].append(hit)
            except: continue
        return grouped

    def _find_amplicons(self, f_hits: List[BlastHit], r_hits: List[BlastHit]) -> List[OffTargetAmplicon]:
        """Strand와 거리 조건을 이용한 증폭 산물 찾기 로직"""
        amplicons = []
        for fh in f_hits:
            for rh in r_hits:
                if fh.sseqid != rh.sseqid: continue
                if fh.strand == rh.strand: continue
                
                valid = False
                start, end = 0, 0
                
                # Forward(+) & Reverse(-) 구도
                if fh.strand == "+" and rh.strand == "-" and fh.genomic_end < rh.genomic_start: 
                    valid, start, end = True, fh.genomic_start, rh.genomic_end
                # Forward(-) & Reverse(+) 구도
                elif fh.strand == "-" and rh.strand == "+" and rh.genomic_end < fh.genomic_start: 
                    valid, start, end = True, rh.genomic_start, fh.genomic_end
                
                if valid:
                    size = end - start
                    if self.criteria.min_amp_size <= size <= self.criteria.max_amp_size:
                        amplicons.append(OffTargetAmplicon(
                            chrom=fh.sseqid, start=start, end=end, 
                            product_size=size, fwd_hit=fh, rev_hit=rh
                        ))
        return amplicons

    def _check_probe_binding(self, cand: OffTargetAmplicon, p_hits: List[BlastHit]) -> bool:
        """증폭 범위 내에 프로브가 결합하는지 확인"""
        for ph in p_hits:
            if ph.sseqid != cand.chrom: continue
            if ph.genomic_start >= cand.start and ph.genomic_end <= cand.end:
                return True
        return False

    def _fetch_sequence_from_genome(self, chrom: str, start: int, end: int) -> str:
        if not self.ref_fasta: return ""
        try:
            return self.ref_fasta.fetch(chrom, start, end).upper()
        except:
            alt = chrom.replace("chr", "") if "chr" in chrom else f"chr{chrom}"
            try: return self.ref_fasta.fetch(alt, start, end).upper()
            except: return ""