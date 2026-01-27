
from typing import List, Dict, Any
from ispcr import calculate_pcr_product, FastaSequence
import pysam

from ..config.schema.qc import QCParams
from .types import OffTargetAmplicon
from ..components import Amplicon

class IsPcrChecker:
    def __init__(self, qc_params: QCParams):
        self.params = qc_params
        # Config에 설정된 Reference Fasta 경로 (hg38.fa 등)
        self.reference_path = qc_params.paths.BLASTN_REF

    def check_amplicon(self, amplicon: Amplicon) -> Dict[str, Any]:
        """
        임시 파일 생성 없이 메모리 상에서 isPCR 수행 (calculate_pcr_product 직접 호출)
        """
        # 1. 필수 라이브러리 체크
        if calculate_pcr_product is None or FastaSequence is None:
            return {"error": "Package 'ispcr' not installed. Run 'pip install ispcr'", "passed": False}
        if pysam is None:
            return {"error": "Package 'pysam' not installed. Run 'pip install pysam'", "passed": False}
        
        # 2. Primer를 FastaSequence 객체로 변환 (메모리 상 생성)
        #    소스 코드상 calculate_pcr_product는 이 객체 타입을 요구함
        fwd_obj = FastaSequence("Forward", amplicon.forward_primer.sequence)
        rev_obj = FastaSequence("Reverse", amplicon.reverse_primer.sequence)

        detected_amplicons = []