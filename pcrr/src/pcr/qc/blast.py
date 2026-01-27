import subprocess
import tempfile
import os
import pysam

from typing import List, Dict, Optional, Tuple, Any
from ispcr import calculate_pcr_product, FastaSequence
from ..config.schema.qc import QCParams
from .types import BlastHit, OffTargetAmplicon
from ..components import Amplicon


class BlastSpecificityChecker:
	OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"

	def __init__(self, qc_params: QCParams):
		self.params = qc_params
		self.criteria = qc_params.get_primer_criteria()
		self.paths = qc_params.paths

	def check_amplicon(self, amplicon: Amplicon) -> Dict[str, Any]:
		"""
		[BLAST + isPCR 하이브리드 QC]
		1. BLAST: 전체 게놈 스캔 -> 잠재적 결합 부위 탐색
		2. isPCR: BLAST Hit 주변부만 잘라내어 정밀 증폭 시뮬레이션
		3. 결과: QC 통과 시 amplicon.is_qc_pass = True 설정
		"""
		fwd = amplicon.forward_primer
		rev = amplicon.reverse_primer
		
		# -------------------------------------------------------
		# 1. BLAST 실행 (Global Search)
		# -------------------------------------------------------
		hits = self._run_blast_pair(
			f_name="Fwd", f_seq=fwd.sequence, 
			r_name="Rev", r_seq=rev.sequence
		)
		f_hits = [h for h in hits if h.qseqid == "Fwd"]
		r_hits = [h for h in hits if h.qseqid == "Rev"]

		# -------------------------------------------------------
		# 2. BLAST 결과로 Amplicon 후보 조합 (Distance 기반)
		# -------------------------------------------------------
		candidates = self._find_amplicons(f_hits, r_hits)

		# -------------------------------------------------------
		# 3. isPCR로 후보 검증 (Local Verification)
		# -------------------------------------------------------
		confirmed_off_targets = []
		intended_target = None

		# isPCR용 프라이머 객체 생성 (반복 생성 방지)
		# 라이브러리가 없을 경우를 대비해 None 처리
		fwd_obj = FastaSequence("Fwd", fwd.sequence) if FastaSequence else None
		rev_obj = FastaSequence("Rev", rev.sequence) if FastaSequence else None

		# Reference 파일 열기 (pysam)
		ref_fasta = None
		try:
			if pysam and os.path.exists(self.paths.BLASTN_REF):
				ref_fasta = pysam.FastaFile(self.paths.BLASTN_REF)
		except Exception as e:
			print(f"⚠️ Warning: Could not open reference genome: {e}")

		# 후보 검증 루프
		for cand in candidates:
			# 3-1. isPCR 검증 (라이브러리가 없으면 BLAST 결과 그대로 신뢰)
			is_real = self._verify_with_ispcr(cand, fwd_obj, rev_obj, ref_fasta)
			
			if is_real:
				if intended_target is None and cand.is_target:
					intended_target = cand
				else:
					confirmed_off_targets.append(cand)
		
		if ref_fasta:
			ref_fasta.close()

		# -------------------------------------------------------
		# 4. 최종 판정 및 Amplicon 객체 업데이트
		# -------------------------------------------------------
		# 통과 조건: 오프타겟이 0개이고, 의도한 타겟은 찾아야 함
		is_passed = (len(confirmed_off_targets) == 0) and (intended_target is not None)
		
		# ✅ 요청하신 기능: 객체에 속성 주입
		amplicon.is_qc_pass = is_passed
		# 추가 정보도 넣어두면 나중에 유용함
		amplicon.off_target_count = len(confirmed_off_targets)

		# -------------------------------------------------------
		# 5. 결과 리턴 (리포트용)
		# -------------------------------------------------------
		alignment_view = None
		if intended_target:
			alignment_view = self._generate_alignment_view(
				amplicon_seq=amplicon.sequence, 
				target_hit=intended_target
			)

		return {
			"passed": is_passed,
			"intended_target_found": intended_target is not None,
			"off_target_count": len(confirmed_off_targets),
			"off_targets": [a.model_dump() for a in confirmed_off_targets],
			"alignment": alignment_view,
			"blast_candidates_count": len(candidates) 
		}

	# -------------------------------------------------------------------------
	# [핵심] BLAST 좌표 기반 isPCR 검증
	# -------------------------------------------------------------------------
	def _verify_with_ispcr(self, cand: OffTargetAmplicon, 
						   fwd_obj: Any, rev_obj: Any, ref_fasta: Any) -> bool:
		"""
		BLAST가 찾은 좌표 주변을 잘라내어 isPCR을 돌림.
		True 리턴 시: 진짜 PCR Product 생성됨.
		"""
		# 라이브러리가 없거나 레퍼런스를 못 열었으면, BLAST 결과를 보수적으로 인정(True)
		if not calculate_pcr_product or not ref_fasta or not fwd_obj:
			return True

		# BLAST 좌표보다 약간 여유있게 가져옴 (Primers가 잘리지 않게 +/- 50bp)
		padding = 50
		# pysam fetch는 0-based indexing
		# cand.start는 1-based (BLAST) -> 0-based 변환 시 -1
		fetch_start = max(0, cand.start - 1 - padding)
		fetch_end = cand.end + padding

		try:
			# 1. Genome에서 해당 조각만 가져옴 (매우 빠름)
			local_seq_str = ref_fasta.fetch(cand.chrom, fetch_start, fetch_end)
			local_seq_obj = FastaSequence(cand.chrom, local_seq_str)

			# 2. 그 조각에 대해서만 isPCR 수행
			result = calculate_pcr_product(
				sequence=local_seq_obj,
				forward_primer=fwd_obj,
				reverse_primer=rev_obj,
				# 조각 내에서의 상대적 길이 조건 (검출만 목적이므로 범위 넓게)
				min_product_length=10, 
				max_product_length=10000,
				header=False,
				cols="all",
				output_file=False
			)
			
			# 결과 문자열이 있으면 증폭 성공
			return bool(result and result.strip())

		except Exception:
			# 에러나면 안전하게 True(위험군)로 간주
			return True

	# -------------------------------------------------------------------------
	# Helper: Alignment View Generator
	# -------------------------------------------------------------------------
	def _generate_alignment_view(self, amplicon_seq: str, target_hit: OffTargetAmplicon) -> Dict[str, str]:
		"""Amplicon과 Genome Reference 비교 시각화"""
		try:
			if not pysam or not os.path.exists(self.paths.BLASTN_REF):
				return {"error": "pysam missing or Ref not found"}

			with pysam.FastaFile(self.paths.BLASTN_REF) as fasta:
				# Genome 서열 추출 (0-based convert)
				ref_seq = fasta.fetch(
					target_hit.chrom, 
					target_hit.start - 1, 
					target_hit.end
				).upper()
			
			amp_seq = amplicon_seq.upper()
			
			match_line = []
			mismatch_count = 0
			
			# 길이 차이 보정 (Visual용)
			min_len = min(len(ref_seq), len(amp_seq))
			
			for i in range(min_len):
				if ref_seq[i] == amp_seq[i]:
					match_line.append("|")
				else:
					match_line.append("*")
					mismatch_count += 1
			
			return {
				"chrom": target_hit.chrom,
				"range": f"{target_hit.start}-{target_hit.end}",
				"mismatch_count": mismatch_count,
				"visualization": (
					f"Ref: {ref_seq[:min_len]}\n"
					f"	 {''.join(match_line)}\n"
					f"Amp: {amp_seq[:min_len]}"
				)
			}
		except Exception as e:
			return {"error": f"Alignment failed: {e}"}

	# -------------------------------------------------------------------------
	# Existing Methods (Run BLAST, Parse, Find Amplicons)
	# -------------------------------------------------------------------------
	def _run_blast_pair(self, f_name: str, f_seq: str, r_name: str, r_seq: str) -> List[BlastHit]:
		"""BLASTN 실행"""
		fasta_content = f">{f_name}\n{f_seq}\n>{r_name}\n{r_seq}\n"
		cmd = [
			self.paths.BLASTN,
			"-task", "blastn-short",
			"-db", self.paths.BLASTN_DB,
			"-outfmt", self.OUTFMT,
			"-num_alignments", str(self.criteria.blast_max_alignments),
			"-num_threads", "2"
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
			if len(cols) < 10: continue

			try:
				hit = BlastHit(
					qseqid=cols[0], sseqid=cols[1],
					pident=float(cols[2]), length=int(cols[3]),
					qstart=int(cols[6]), qend=int(cols[7]),
					sstart=int(cols[8]), send=int(cols[9])
				)
				if (hit.pident >= self.criteria.blast_identity_threshold and 
					hit.length >= self.criteria.blast_length_threshold):
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
				
				valid_orientation = False
				if fh.strand == "+" and rh.strand == "-" and fh.genomic_end < rh.genomic_start:
					valid_orientation = True
				elif fh.strand == "-" and rh.strand == "+":
					# Circular DNA or unusual case, usually R end < F start for standard PCR
					# But for simple linear logic:
					if rh.genomic_end < fh.genomic_start: valid_orientation = True
				
				if not valid_orientation: continue

				start = min(fh.genomic_start, rh.genomic_start)
				end = max(fh.genomic_end, rh.genomic_end)
				product_size = end - start + 1
				
				# BLAST 단계에서는 아주 넓은 범위로 잡아두는 것이 좋음 (isPCR에서 거르면 됨)
				# 여기서는 Config 기준을 따름
				if self.criteria.min_amp_bp <= product_size <= self.criteria.max_amp_bp:
					amplicons.append(
						OffTargetAmplicon(
							chrom=fh.sseqid,
							start=start,
							end=end,
							product_size=product_size,
							fwd_hit=fh,
							rev_hit=rh
						)
					)
		return amplicons