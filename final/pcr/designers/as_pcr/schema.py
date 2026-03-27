from typing import Dict, Optional, List, Any
from pydantic import BaseModel, Field, model_validator

from pcr.components.amplicon import Amplicon
from pcr.designers.base.schema import BaseDesignInput, BaseDesignOutput

class ASPCRAmplicon(Amplicon):
    """
    AS-PCR 전용 앰플리콘 모델. 
    일반 Amplicon에 AS-PCR 분석 및 시각화용 메타데이터 필드를 추가합니다.
    """
    allele_type: Optional[str] = None
    set_id: Optional[str] = None
    fixed_prime: Optional[str] = None
    alignment_visual: Optional[List[str]] = None

class ASPCRDesignInput(BaseDesignInput):
    """
    AS-PCR 전용 입력 객체.
    반복되는 파라미터와 변환 로직은 BaseDesignInput을 상속받아 자동으로 처리하고,
    AS-PCR 특화 설정만 추가로 주입합니다.
    """
    design_name: str = "AS-PCR_Design"
    sequence: str = "" 
    reference_genome: str = "hg38"
    top_k: int = 5
    
    template_genomic_start: int = 0
    template_genomic_end: int = 0

    # 🔥 AS-PCR 특화 설정 (자식 클래스의 고유 속성)
    fixed_prime: str = "forward"
    mismatch_pos: int = 3
    mismatch_intensity: str = "strong"

    @model_validator(mode='before')
    @classmethod
    def parse_brackets_and_validate(cls, data: Any) -> Any:
        """
        API에서 넘어온 데이터를 Base 파서가 이해할 수 있는 규격으로 이름만 맞춰줍니다.
        (대괄호 파싱 등은 이후 Base의 부모 Validator가 이어서 처리하거나 무시합니다)
        """
        if isinstance(data, dict):
            data['name'] = data.get('design_name', 'AS-PCR_Design')
            data['reference_name'] = data.get('reference_genome', 'hg38')
            
            if 'sequence' in data:
                data['template_sequence'] = data['sequence']
                
            data['target_start'] = data.get('target_start', 0)
            data['target_end'] = data.get('target_end', 0)
            
        return data

    def to_core_qc_overrides(self) -> Dict[str, Any]:
        """
        [Adapter] 
        부모의 공통 QC 변환 로직을 그대로 호출한 뒤, AS-PCR 전용 블록만 덧붙입니다.
        """
        overrides = super().to_core_qc_overrides()
        
        # Base에서 만들어준 딕셔너리에 AS-PCR 고유 설정 결합
        overrides["as_pcr"] = {
            "fixed_prime": self.fixed_prime,
            "mismatch_pos": self.mismatch_pos,
            "mismatch_intensity": self.mismatch_intensity
        }
        
        return overrides

class ASPCRDesignOutput(BaseDesignOutput):
    """
    AS-PCR 전용 출력 스키마.
    향후 세트(Set) 단위의 그룹화 정보나 AS-PCR 특화 메타데이터를 담을 때 확장합니다.
    """
    pass