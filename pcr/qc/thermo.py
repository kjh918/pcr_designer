"""
pcr/qc/thermo.py
순수 열역학(Secondary Structure) 검증 도구.
AmpliconQCStatus의 add_result 인터페이스를 사용하여 결과를 누적합니다.
"""
import primer3
from typing import Dict, Any, List, Union
from pcr.config.schema.qc import QCCriteria, AmpliconQCStatus
from pcr.components.amplicon import Amplicon
from pcr.components.primer import Primer, Probe

class ThermoChecker:
    def __init__(self, criteria: QCCriteria):
        self.criteria = criteria

    def run(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        """
        [Standard Interface] 리스트 전체를 검사하여 결과를 각 앰플리콘의 qc_status에 누적합니다.
        """
        for amp in amplicons:
            # 1. 앰플리콘에 qc_status 객체가 없으면 초기화
            if not hasattr(amp, "qc_status") or amp.qc_status is None:
                amp.qc_status = AmpliconQCStatus()

            # 2. 단일 검사 수행
            result = self._check_single(amp)
            
            # 3. 🔥 최신 Workflow 구조(add_result) 반영
            # fail_reason 문자열을 리스트로 변환하여 전달
            messages = [m.strip() for m in result["fail_reason"].split(",") if m.strip()]
            
            amp.qc_status.add_result(
                module_name="thermo",
                is_pass=result["passed"],
                messages=messages,
                metrics=result["data"]
            )
            
            # 4. 하위 호환성을 위한 속성 업데이트 (필요 시)
            amp.is_qc_pass = amp.qc_status.is_pass
            # 모든 모듈의 fail_reasons를 합쳐서 qc_log 업데이트
            amp.qc_log = " | ".join(amp.qc_status.fail_reasons)

        return amplicons

    def check(self, amplicon: Amplicon) -> Dict[str, Any]:
        """
        [Bridge Interface] 단일 앰플리콘의 상세 검사 데이터를 반환합니다.
        """
        return self._check_single(amplicon)

    def _check_single(self, amplicon: Amplicon) -> Dict[str, Any]:
        """정교한 열역학 검사 로직 (Hairpin, Homo/Hetero Dimer)"""
        details = {}
        all_pass = True
        fail_reasons = []

        # 1. Individual Oligo Check (Hairpin & Homodimer)
        f_res = self._check_oligo(amplicon.forward, "fwd")
        details.update(f_res)
        if not f_res["fwd_pass"]:
            all_pass = False
            fail_reasons.append(f"Fwd_Thermo({f_res.get('fwd_fail_reason')})")

        r_res = self._check_oligo(amplicon.reverse, "rev")
        details.update(r_res)
        if not r_res["rev_pass"]:
            all_pass = False
            fail_reasons.append(f"Rev_Thermo({r_res.get('rev_fail_reason')})")

        if amplicon.probe:
            p_res = self._check_oligo(amplicon.probe, "probe")
            details.update(p_res)
            if not p_res["probe_pass"]:
                all_pass = False
                fail_reasons.append(f"Probe_Thermo({p_res.get('probe_fail_reason')})")
        else:
            details.update({
                "probe_hairpin_dg": 0.0, 
                "probe_homodimer_dg": 0.0, 
                "probe_pass": True
            })

        # 2. Heterodimer Check (Interactions)
        fr_res = self._check_hetero(amplicon.forward.sequence, amplicon.reverse.sequence, "hetero_fr")
        details.update(fr_res)
        if not fr_res["hetero_fr_pass"]:
            all_pass = False
            fail_reasons.append("Hetero_FR_Dimer")

        if amplicon.probe:
            fp_res = self._check_hetero(amplicon.forward.sequence, amplicon.probe.sequence, "hetero_fp")
            rp_res = self._check_hetero(amplicon.reverse.sequence, amplicon.probe.sequence, "hetero_rp")
            details.update(fp_res); details.update(rp_res)
            if not fp_res["hetero_fp_pass"]: all_pass = False; fail_reasons.append("Hetero_FP_Dimer")
            if not rp_res["hetero_rp_pass"]: all_pass = False; fail_reasons.append("Hetero_RP_Dimer")
        else:
            details.update({
                "hetero_fp_dg": 0.0, "hetero_fp_pass": True, 
                "hetero_rp_dg": 0.0, "hetero_rp_pass": True
            })

        return {
            "passed": all_pass,
            "fail_reason": ", ".join(fail_reasons) if fail_reasons else "",
            "data": details
        }

    def _check_oligo(self, oligo: Union[Primer, Probe], prefix: str) -> Dict[str, Any]:
        """개별 올리고의 Hairpin 및 Homodimer 계산"""
        seq = oligo.sequence
        tm = oligo.tm
        reasons = []
        # QCCriteria 스키마의 값을 사용 (pcr/config/schema/qc.py 참조)
        min_hp = self.criteria.hairpin_min_dg
        min_hd = self.criteria.homodimer_min_dg

        try:
            hp = primer3.calc_hairpin(seq, temp_c=tm)
            hp_dg = hp.dg / 1000.0 if hp.structure_found else 0.0
        except Exception: hp_dg = 0.0

        try:
            hd = primer3.calc_homodimer(seq, temp_c=tm)
            hd_dg = hd.dg / 1000.0 if hd.structure_found else 0.0
        except Exception: hd_dg = 0.0

        hp_pass = (hp_dg >= min_hp)
        hd_pass = (hd_dg >= min_hd)
        if not hp_pass: reasons.append("Hairpin")
        if not hd_pass: reasons.append("Homodimer")

        return {
            f"{prefix}_hairpin_dg": round(hp_dg, 2),
            f"{prefix}_homodimer_dg": round(hd_dg, 2),
            f"{prefix}_pass": hp_pass and hd_pass,
            f"{prefix}_fail_reason": "/".join(reasons)
        }

    def _check_hetero(self, s1: str, s2: str, label: str) -> Dict[str, Any]:
        """두 올리고 간의 Heterodimer 계산"""
        min_het = self.criteria.heterodimer_min_dg
        try:
            het = primer3.calc_heterodimer(s1, s2)
            het_dg = het.dg / 1000.0 if het.structure_found else 0.0
        except Exception: het_dg = 0.0
        return {
            f"{label}_dg": round(het_dg, 2), 
            f"{label}_pass": (het_dg >= min_het)
        }