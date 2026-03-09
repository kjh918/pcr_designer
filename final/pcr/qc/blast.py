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
    
    OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
    
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
                amp.is_qc_pass = result_data["passed"]
                
                if amp.is_qc_pass:
                    amp.qc_log = "Passed: Specific target confirmed."
                else:
                    # 실패 사유를 상세히 기록 (프론트에서 보여줄 용도)
                    amp.qc_log = result_data.get("reason", "Failed: Specificity issue.")

        except Exception as e:
            print(f"❌ BLAST Execution Failed: {e}")
            for amp in amplicons:
                amp.is_qc_pass = False
                amp.qc_log = f"Error: BLAST analysis failed ({str(e)})"

        # 🎯 필터링 없이 전체 리스트 반환 (프론트에서 is_qc_pass로 구분해서 그리게 함)
        return amplicons
    def _create_binding_report(self, cand: OffTargetAmplicon, amplicon: Amplicon, p_hits: List[BlastHit], is_target: bool = False) -> Dict[str, Any]:
        """
        [Core] 타겟/오프타겟 영역의 Fwd/Rev/Probe 바인딩 상세 정보를 시각화 데이터로 구성합니다.
        (객체를 직접 변조하지 않고 is_target 파라미터를 받아 안전하게 처리합니다)
        """
        
        # 1. Forward Primer Alignment
        f_report = self._generate_alignment_block(
            name="Forward Primer", query_seq=cand.fwd_hit.qseq, subject_seq=cand.fwd_hit.sseq, pident=cand.fwd_hit.pident
        )

        # 2. Reverse Primer Alignment
        r_report = self._generate_alignment_block(
            name="Reverse Primer", query_seq=cand.rev_hit.qseq, subject_seq=cand.rev_hit.sseq, pident=cand.rev_hit.pident
        )

        # 3. Probe Alignment
        p_report = None
        target_p_hit = next((ph for ph in p_hits if ph.genomic_start >= cand.start and ph.genomic_end <= cand.end), None)
        if target_p_hit:
            p_report = self._generate_alignment_block(
                name="TaqMan Probe", query_seq=target_p_hit.qseq, subject_seq=target_p_hit.sseq, pident=target_p_hit.pident
            )

        return {
            "location": f"{cand.chrom}:{cand.start}-{cand.end}",
            "product_size": cand.product_size,
            "is_target": is_target, # 🔥 매개변수로 받아 안전하게 주입
            "fwd": f_report,
            "rev": r_report,
            "probe": p_report,
            "full_sequence": getattr(cand, "reference_sequence", "") 
        }

    def _evaluate_amplicon(self, amplicon: Amplicon, f_hits: List[BlastHit], r_hits: List[BlastHit], p_hits: List[BlastHit]) -> Dict[str, Any]:
        """[Sequence-Only Mode] 객체 변조 에러를 해결하고 프로브 유무를 동적으로 처리하는 정밀 판정 로직"""
        candidates = self._find_amplicons(f_hits, r_hits)
        
        target_signals = []      
        off_target_signals = []  
        amplification_only = []  

        f_len = len(amplicon.forward.sequence)
        r_len = len(amplicon.reverse.sequence)
        
        # 🔥 추가: 현재 앰플리콘에 프로브가 포함되어 있는지 확인 (AS-PCR은 프로브 없음!)
        has_probe = amplicon.probe is not None

        min_cov = getattr(self.criteria, "min_query_coverage", 0.8)
        min_id = getattr(self.criteria, "min_identity_threshold", 90.0)
        tolerance = getattr(self.criteria, "end_match_tolerance", 1)

        for idx, cand in enumerate(candidates):
            f_cov = cand.fwd_hit.length / f_len
            f_3end = (cand.fwd_hit.qend >= f_len - tolerance)
            f_valid = (f_cov >= min_cov) and f_3end and (cand.fwd_hit.pident >= min_id)

            r_cov = cand.rev_hit.length / r_len
            r_3end = (cand.rev_hit.qstart <= 1 + tolerance)
            r_valid = (r_cov >= min_cov) and r_3end and (cand.rev_hit.pident >= min_id)

            if not (f_valid and r_valid):
                continue

            # 🔥 프로브가 없는 Assay(AS-PCR)면 묻지도 따지지도 않고 결합(True)한 것으로 간주!
            probe_binds = self._check_probe_binding(cand, p_hits) if has_probe else True
            
            # 객체 조작 싹 지우고, is_target 변수로 내려보냄
            is_target = (idx == 0)
            binding_report = self._create_binding_report(cand, amplicon, p_hits if has_probe and probe_binds else [], is_target)

            if is_target:
                if probe_binds:
                    target_signals.append(binding_report)
            else:
                if probe_binds:
                    off_target_signals.append(binding_report)
                else:
                    amplification_only.append(binding_report)

        # 최종 QC 판정 (타겟 1개 필수, 오프타겟 및 노이즈 0개 필수)
        is_passed = (len(target_signals) == 1) and \
                    (len(off_target_signals) == 0) and \
                    (len(amplification_only) == 0)
            
        # 프론트에 보여줄 에러 사유 로직 분리
        reason = ""
        if not is_passed:
            if len(target_signals) == 0:
                reason = "Failed: Main target amplification failed" + (" or Probe missing." if has_probe else ".")
            elif len(off_target_signals) > 0 or len(amplification_only) > 0:
                reason = f"Failed: Off-target detected ({len(off_target_signals)} probe-binds, {len(amplification_only)} amp-only)."

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
    
    def _generate_alignment_block(self, name: str, query_seq: str, subject_seq: str, pident: float) -> Dict[str, str]:
        """미스매치 시각화를 위한 문자열 블록 생성"""
        match_line = ""
        for q, s in zip(query_seq, subject_seq):
            if q == s:
                match_line += "|"
            elif q == "-" or s == "-": # 인델(Gap) 처리
                match_line += " "
            else: # 미스매치
                match_line += "." # 또는 공백
        
        return {
            "label": name,
            "query": query_seq,
            "match": match_line,
            "subject": subject_seq,
            "identity": f"{pident}%"
        }
    
    def _evaluate_amplicon(self, amplicon: Amplicon, f_hits: List[BlastHit], r_hits: List[BlastHit], p_hits: List[BlastHit]) -> Dict[str, Any]:
        candidates = self._find_amplicons(f_hits, r_hits)
        
        target_signals = []      
        off_target_signals = []  
        amplification_only = []  

        f_len = len(amplicon.forward.sequence)
        r_len = len(amplicon.reverse.sequence)
        
        # 🔥 추가: 현재 앰플리콘에 프로브가 포함되어 있는지 확인 (AS-PCR은 프로브 없음!)
        has_probe = amplicon.probe is not None

        min_cov = getattr(self.criteria, "min_query_coverage", 0.8)
        min_id = getattr(self.criteria, "min_identity_threshold", 90.0)
        tolerance = getattr(self.criteria, "end_match_tolerance", 1)

        for idx, cand in enumerate(candidates):
            f_cov = cand.fwd_hit.length / f_len
            f_3end = (cand.fwd_hit.qend >= f_len - tolerance)
            f_valid = (f_cov >= min_cov) and f_3end and (cand.fwd_hit.pident >= min_id)

            r_cov = cand.rev_hit.length / r_len
            r_3end = (cand.rev_hit.qstart <= 1 + tolerance)
            r_valid = (r_cov >= min_cov) and r_3end and (cand.rev_hit.pident >= min_id)

            if not (f_valid and r_valid):
                continue

            # 🔥 프로브가 없는 Assay(AS-PCR)면 묻지도 따지지도 않고 결합(True)한 것으로 간주!
            probe_binds = self._check_probe_binding(cand, p_hits) if has_probe else True
            
            # 객체 조작 싹 지우고, is_target 변수로 내려보냄
            is_target = (idx == 0)
            binding_report = self._create_binding_report(cand, amplicon, p_hits if has_probe and probe_binds else [], is_target)

            if is_target:
                if probe_binds:
                    target_signals.append(binding_report)
            else:
                if probe_binds:
                    off_target_signals.append(binding_report)
                else:
                    amplification_only.append(binding_report)

        # 최종 QC 판정 (타겟 1개 필수, 오프타겟 및 노이즈 0개 필수)
        is_passed = (len(target_signals) == 1) and \
                    (len(off_target_signals) == 0) and \
                    (len(amplification_only) == 0)
            
        print(len(target_signals), len(off_target_signals), len(amplification_only))
        # 프론트에 보여줄 에러 사유 로직 분리
        reason = ""
        if not is_passed:
            if len(target_signals) == 0:
                reason = "Failed: Main target amplification failed" + (" or Probe missing." if has_probe else ".")
            elif len(off_target_signals) > 0 or len(amplification_only) > 0:
                reason = f"Failed: Off-target detected ({len(off_target_signals)} probe-binds, {len(amplification_only)} amp-only)."

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
    
    def _is_potential_amplification(self, hit: BlastHit, primer_len: int) -> bool:
        """
        YAML 설정값을 기반으로 실제 증폭 가능성을 판정합니다.
        """
        # 1. YAML에서 값 로드 (기본값 설정으로 안전성 확보)
        min_cov = getattr(self.criteria, "min_query_coverage", 0.8)
        min_id = getattr(self.criteria, "min_identity_threshold", 90.0)
        tolerance = getattr(self.criteria, "end_match_tolerance", 1)

        # 2. Query Coverage 체크
        coverage = hit.length / primer_len
        if coverage < min_cov:
            return False

        # 3. 3' 말단(3'-end) 결합 체크
        is_3end_match = False
        if hit.qseqid == "Fwd":
            # Forward 프라이머의 끝(qend)이 프라이머 전체 길이 근처인지
            is_3end_match = (hit.qend >= primer_len - tolerance)
        elif hit.qseqid == "Rev":
            # Reverse 프라이머는 qstart가 시작점(1) 근처여야 3' 말단 결합
            is_3end_match = (hit.qstart <= 1 + tolerance)

        # 4. 종합 판정
        # Identity가 설정값(90%) 이상이고 3' 말단이 붙어있어야 '진짜 증폭 위험'
        return is_3end_match and (hit.pident >= min_id)
    
    def _parse_and_group_hits(self, stdout: str) -> Dict[str, List[BlastHit]]:
        grouped = {}
        for line in stdout.strip().splitlines():
            cols = line.split("\t")
            # 🚨 qseq, sseq를 포함하려면 최소 14개의 컬럼이 필요합니다.
            if len(cols) < 14: 
                continue
            
            raw_qseqid = cols[0]
            if "|" not in raw_qseqid: 
                continue
            amp_id, oligo_type = raw_qseqid.split("|", 1)

            try:
                sstart, send = int(cols[8]), int(cols[9])
                
                # 1. BlastHit 객체 생성 (서열 정보 포함)
                hit = BlastHit(
                    qseqid=oligo_type, 
                    sseqid=cols[1], 
                    pident=float(cols[2]), 
                    length=int(cols[3]), 
                    qstart=int(cols[6]), 
                    qend=int(cols[7]),
                    sstart=sstart, 
                    send=send,
                    # 🔥 추가: BLAST가 매칭한 Query 및 Subject 서열
                    qseq=cols[12], 
                    sseq=cols[13]
                )
                # 2. 필터링 및 그룹화
                # 팁: 프로브(Probe)는 프라이머보다 더 낮은 identity에서도 
                # 바인딩할 수 있으므로, oligo_type에 따라 기준을 다르게 줄 수도 있습니다.
                min_id = self.criteria.min_identity
                if oligo_type == "Probe":
                    # 프로브용 별도 기준이 있다면 적용 (예: 80.0)
                    min_id = getattr(self.criteria, "probe_min_identity", min_id)

                if hit.pident >= min_id:
                    if amp_id not in grouped: 
                        grouped[amp_id] = []
                    grouped[amp_id].append(hit)
                    
            except Exception as e:
                # 로깅을 추가하면 파싱 에러 발생 시 원인 파악이 쉽습니다.
                # print(f"Error parsing line: {e}")
                continue
                
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