from typing import Optional, Dict, Any
from pydantic import BaseModel, Field, model_validator
from pcr.designers.base.schema import BaseDesignInput, BaseDesignOutput

class QPCRDesignInput(BaseDesignInput):
    """
    qPCR 전용 입력 객체. 
    API 통신부터 Factory 주입까지 이 객체 하나로 처리합니다.
    """
    # 프론트엔드와 매칭되는 필드들
    design_name: str = "qPCR_Design"
    sequence: str = "" 
    reference_genome: str = "hg38"
    top_k: int = 10
    
    # 🔥 핵심 버그 수정: TypeError ('>' not supported between 'NoneType' and 'int') 방지
    # getattr()에서 None이 반환되어 designer.py에서 에러가 터지는 것을 막기 위해 기본값을 0으로 강제합니다.
    template_genomic_start: int = 0
    template_genomic_end: int = 0
    
    # 계층형 파라미터 (app.js 페이로드 호환성 유지)
    pcr_params: Dict[str, Any] = Field(default_factory=dict)
    qc_criteria: Dict[str, Any] = Field(default_factory=dict)
    
    # 이전 버전의 payload 지원용
    amplicon: Dict[str, Any] = Field(default_factory=dict)
    primer: Dict[str, Any] = Field(default_factory=dict)
    probe: Dict[str, Any] = Field(default_factory=dict)
    
    # qPCR 특화 설정
    require_probe: bool = True
    target_dye: Optional[str] = "FAM"
    quencher: Optional[str] = "BHQ1"

    @model_validator(mode='before')
    @classmethod
    def parse_brackets_and_validate(cls, data: Any) -> Any:
        """
        🔥 부모(BaseDesignInput)의 단순 파싱 로직을 덮어쓰고(Override) 무효화시킵니다!
        프론트엔드에서 날아오는 API JSON(sequence)을 파싱하여 [REF, ALT] 형식을 유지하고, 
        코어 엔진이 요구하는 필수 필드(template_sequence, target_start 등)를 안전하게 충족시킵니다.
        """
        if isinstance(data, dict):
            # 1. API 파라미터 매핑 (이름 통일)
            if 'design_name' in data:
                data['name'] = data['design_name']
            if 'reference_genome' in data:
                data['reference_name'] = data['reference_genome']

            # 2. 서열 파싱 (sequence -> template_sequence / reference_sequence)
            seq_input = data.get('sequence', '') or data.get('template_sequence', '')
            
            # 괄호가 있는 경우 [REF, ALT] 정밀 파싱 수행
            if seq_input and '[' in seq_input and ']' in seq_input:
                seq = seq_input.replace("\n", "").upper()
                template_seq = ""
                reference_seq = ""
                target_indices = []
                
                in_target = False
                bracket_content = ""
                current_idx = 0
                
                for char in seq:
                    if char == '[':
                        in_target = True
                        bracket_content = ""
                    elif char == ']':
                        in_target = False
                        
                        # 콤마(,) 또는 슬래시(/)를 기준으로 분리
                        if ',' in bracket_content or '/' in bracket_content:
                            delimiter = ',' if ',' in bracket_content else '/'
                            parts = [p.replace(" ", "") for p in bracket_content.split(delimiter)]
                            ref_allele = parts[0]
                            alt_allele = parts[1] if len(parts) > 1 else parts[0]
                        else:
                            ref_allele = bracket_content.replace(" ", "")
                            alt_allele = bracket_content.replace(" ", "")
                        
                        for _ in range(len(alt_allele)):
                            target_indices.append(current_idx)
                            current_idx += 1
                        
                        template_seq += alt_allele
                        reference_seq += ref_allele
                    else:
                        if in_target:
                            bracket_content += char
                        else:
                            if char != " ":
                                template_seq += char
                                reference_seq += char
                                current_idx += 1

                # 부모 클래스가 요구하는 필수값 주입
                data['template_sequence'] = template_seq
                data['reference_sequence'] = reference_seq
                data['sequence'] = seq_input 
                
                if target_indices:
                    data['target_start'] = min(target_indices)
                    data['target_end'] = max(target_indices)
                else:
                    raise ValueError("유효한 타겟 좌표를 찾을 수 없습니다.")
            
            # 괄호가 없는 경우 좌표 필수 확인 (팩토리 내부 통신용)
            else:
                if data.get('target_start') is None or data.get('target_end') is None:
                    raise ValueError("서열에 대괄호('[', ']')로 타겟을 지정하거나 좌표를 명시해야 합니다.")
                if data['target_start'] > data['target_end']:
                    raise ValueError(f"target_start({data['target_start']})는 target_end({data['target_end']})보다 클 수 없습니다.")

        return data

    def to_core_pcr_params(self) -> Dict[str, Any]:
        """[Adapter] 프론트 JSON 딕셔너리 -> 코어 엔진이 이해하는 PCR 설계 파라미터 변환"""
        # 최신 app.js 형식을 우선 지원
        if self.pcr_params:
            return self.pcr_params
            
        # 예전 app.js(amplicon, primer, probe 분리) 형식도 안전하게 지원
        return {
            "primer_kwargs": {
                "n_candidates": 3,
                "min_amplicon_length": self.amplicon.get("min_length", 60),
                "max_amplicon_length": self.amplicon.get("max_length", 150),
                "min_length": self.primer.get("min_length", 20),
                "opt_length": self.primer.get("opt_length", 25),
                "max_length": self.primer.get("max_length", 30),
                "min_tm": self.primer.get("min_tm", 55.0),
                "opt_tm": self.primer.get("opt_tm", 60.0),
                "max_tm": self.primer.get("max_tm", 65.0),
                "min_gc": self.primer.get("min_gc", 35.0),
                "opt_gc": self.primer.get("opt_gc", 50.0),
                "max_gc": self.primer.get("max_gc", 65.0),
            },
            "probe_kwargs": {
                "n_candidates": 5,
                "min_length": self.probe.get("min_length", 20),
                "opt_length": self.probe.get("opt_length", 25),
                "max_length": self.probe.get("max_length", 30),
                "min_tm": self.probe.get("min_tm", 65.0),
                "opt_tm": self.probe.get("opt_tm", 67.0),
                "max_tm": self.probe.get("max_tm", 70.0),
                "min_gc": self.probe.get("min_gc", 35.0),
                "opt_gc": self.probe.get("opt_gc", 50.0),
                "max_gc": self.probe.get("max_gc", 65.0),
                "max_probe_poly_g": self.probe.get("max_poly_g", 3),
                "max_probe_3_end_gc": self.probe.get("max_3_end_gc", 2),
                "avoid_5_prime_g": self.probe.get("avoid_5_prime_g", True)
            }
        }

    def to_core_qc_overrides(self) -> Dict[str, Any]:
        """[Adapter] 프론트 JSON 딕셔너리 -> 코어 엔진이 이해하는 QC 설정 변환"""
        if self.qc_criteria and "hairpin_min_dg" in self.qc_criteria:
            return self.qc_criteria

        thermo = self.qc_criteria.get("thermodynamics", {})
        blast = self.qc_criteria.get("blast", {})
        amp = self.qc_criteria.get("amplicon", {})
        oligo = self.qc_criteria.get("oligo", {})
        
        return {
            "hairpin_min_dg": thermo.get("hairpin_min_dg", -5.0),
            "homodimer_min_dg": thermo.get("homodimer_min_dg", -6.0),
            "heterodimer_min_dg": thermo.get("heterodimer_min_dg", -6.0),
            "min_identity": blast.get("min_identity", 90.0),
            "min_identity_threshold": blast.get("min_identity", 90.0),
            "min_hit_length": blast.get("min_hit_length", 13),
            "blast_max_alignments": blast.get("max_alignments", 50),
            "min_amp_size": amp.get("min_size", 50),
            "max_amp_size": amp.get("max_size", 300),
            "use_ispcr_check": amp.get("use_ispcr", False),
            "primer": {"max_diff_tm": oligo.get("max_tm_diff", 3.0)}
        }

class QPCRDesignOutput(BaseDesignOutput):
    probe_designed_count: int = 0