from typing import List, Dict, Any, Optional
import os

# 외부 라이브러리 예외처리
try:
    from ispcr import calculate_pcr_product, FastaSequence
    import pysam
except ImportError:
    calculate_pcr_product = None
    FastaSequence = None
    pysam = None

# ✅ 수정 1: QCParams 대신 AppConfig 임포트
from ..config.schema.root import AppConfig
from ..components.amplicon import Amplicon

class IsPcrChecker:
    """
    In-Silico PCR 검증기 (Local Verification)
    Reference Genome의 특정 구간(Target 주변)을 추출하여
    실제 PCR 증폭이 일어나는지, 사이즈는 정확한지 검증합니다.
    """
    
    # ✅ 수정 2: 초기화 시 config: AppConfig 받기
    def __init__(self, config: AppConfig):
        self.config = config
        
        # ✅ 수정 3: 변경된 경로 구조에 맞게 접근 (qc_tools.paths.blast_ref_path)
        # schema/qc.py의 QCPaths 정의를 따름
        self.reference_path = config.qc_tools.paths.blast_ref_path

    def check_amplicon(self, amplicon: Amplicon) -> Dict[str, Any]:
        """
        메모리 상에서 타겟 구간에 대한 isPCR 수행
        """
        # 1. 필수 라이브러리 및 파일 체크
        if calculate_pcr_product is None or FastaSequence is None:
            return {"error": "Package 'ispcr' missing", "passed": False}
        if pysam is None:
            return {"error": "Package 'pysam' missing", "passed": False}
        
        # 경로 유효성 체크
        if not self.reference_path or not os.path.exists(self.reference_path):
            return {"error": f"Reference file not found: {self.reference_path}", "passed": False}
        
        # 2. Primer 객체 생성 (ispcr 라이브러리용)
        f_seq = amplicon.forward.sequence
        r_seq = amplicon.reverse.sequence
        
        fwd_obj = FastaSequence("Forward", f_seq)
        rev_obj = FastaSequence("Reverse", r_seq)
        detected_products = []
        
        # 3. Pysam을 이용한 Reference 구간 추출 (Local Fetching)
        try:
            with pysam.FastaFile(self.reference_path) as ref_fasta:
                # Chromosome 이름 (hg38 등)
                chrom = amplicon.reference_id 
                
                # 검색 범위 설정: 타겟 좌표 앞뒤로 여유(Padding)를 둠
                padding = 100
                fetch_start = max(0, amplicon.target_start_index - padding)
                fetch_end = amplicon.target_end_index + padding
                
                # 해당 구간 서열 가져오기
                template_seq_str = ref_fasta.fetch(chrom, fetch_start, fetch_end)
                
                # ispcr용 템플릿 객체 생성
                # 이름 포맷: chrom:start-end
                template_obj = FastaSequence(f"{chrom}:{fetch_start}-{fetch_end}", template_seq_str)
                
                # 4. calculate_pcr_product 실행
                result = calculate_pcr_product(
                    sequence=template_obj,
                    forward_primer=fwd_obj,
                    reverse_primer=rev_obj,
                    min_product_length=30,
                    max_product_length=500,
                    header=False,
                    cols="all",
                    output_file=False
                )
                print(result)
                exit()
                if result and result.strip():
                    # 증폭 성공 시 정보 저장
                    product_seq = result.strip()
                    detected_products.append({
                        "chrom": chrom,
                        "size": len(product_seq),
                        "sequence": product_seq,
                        "is_target": True
                    })

        except Exception as e:
            # pysam fetch 실패 등
            return {"error": f"isPCR execution failed: {str(e)}", "passed": False}

        # 5. 결과 판정
        is_passed = len(detected_products) > 0
        
        return {
            "passed": is_passed,
            "detected_products": detected_products,
            "amplicon_size_check": detected_products[0]['size'] if is_passed else 0,
            "msg": "Target verified" if is_passed else "No amplification at target region"
        }