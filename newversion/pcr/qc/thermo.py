import primer3
from typing import Dict, Any
from ..config.schema.qc import QCCriteria, PrimerQCCriteria, ProbeQCCriteria
from ..components.amplicon import Amplicon, Primer, Probe

class ThermoChecker:
    """
    물리적 성질(Hairpin, Dimer dG)을 검사하는 클래스.
    Primer3-py의 계산 기능을 사용하여 dG(Gibbs Free Energy)를 산출합니다.
    """
    def __init__(self, criteria: QCCriteria):
        # QCCriteria 루트 객체를 저장 (여기에 hairpin_min_dg 등이 있음)
        self.criteria = criteria
        # 하위 기준도 편의상 저장
        self.primer_criteria = criteria.primer
        self.probe_criteria = criteria.probe

    def check(self, amplicon: Amplicon) -> Dict[str, Any]:
        """
        Amplicon 내의 모든 올리고(F, R, P)에 대해 물리적 안정성을 검사합니다.
        """
        details = {}
        all_pass = True
        fail_reasons = []

        # 1. Forward Primer Check
        # ✅ 수정: self.primer_criteria가 아니라 self.criteria(Root)를 넘김
        f_res = self._check_oligo(amplicon.forward, self.criteria, "fwd")
        details.update(f_res)
        if not f_res["fwd_pass"]: 
            all_pass = False
            fail_reasons.append(f"Fwd_Thermo({f_res.get('fwd_fail_reason', '')})")

        # 2. Reverse Primer Check
        r_res = self._check_oligo(amplicon.reverse, self.criteria, "rev")
        details.update(r_res)
        if not r_res["rev_pass"]: 
            all_pass = False
            fail_reasons.append(f"Rev_Thermo({r_res.get('rev_fail_reason', '')})")

        # 3. Primer Heterodimer (Forward + Reverse)
        # ✅ 수정: self.criteria(Root) 넘김 (heterodimer_min_dg 접근용)
        fr_res = self._check_hetero(
            amplicon.forward.sequence, 
            amplicon.reverse.sequence, 
            self.criteria, 
            "hetero_fr"
        )
        details.update(fr_res)
        if not fr_res["hetero_fr_pass"]: 
            all_pass = False
            fail_reasons.append("Hetero_FR")

        # 4. Probe Check (Optional)
        if amplicon.probe:
            p_res = self._check_oligo(amplicon.probe, self.criteria, "probe")
            details.update(p_res)
            if not p_res["probe_pass"]: 
                all_pass = False
                fail_reasons.append(f"Probe_Thermo({p_res.get('probe_fail_reason', '')})")

            # 4-1. Probe Heterodimer (F+P, R+P)
            fp_res = self._check_hetero(amplicon.forward.sequence, amplicon.probe.sequence, self.criteria, "hetero_fp")
            rp_res = self._check_hetero(amplicon.reverse.sequence, amplicon.probe.sequence, self.criteria, "hetero_rp")
            
            if not (fp_res["hetero_fp_pass"] and rp_res["hetero_rp_pass"]):
                all_pass = False
                fail_reasons.append("Hetero_Probe")
            
            details.update(fp_res)
            details.update(rp_res)

        # 결과 통합
        return {
            "passed": all_pass,
            "reason": ", ".join(fail_reasons) if fail_reasons else "Pass",
            "data": details
        }

    def _check_oligo(self, oligo: Primer, criteria: QCCriteria, prefix: str) -> Dict[str, Any]:
        """개별 올리고뉴클레오티드의 Hairpin 및 Homodimer 검사"""
        seq = oligo.sequence
        reasons = []
        
        # ✅ 수정: snake_case 함수명 사용 (Warning 해결)
        # calc_hairpin은 결과를 객체로 반환하며, .dg 속성을 가짐
        try:
            hp = primer3.calc_hairpin(seq, temp_c=oligo.tm)
            hp_dg = hp.dg / 1000.0 if hp.structure_found else 0.0
        except:
            hp_dg = 0.0
        
        # ✅ 수정: snake_case 함수명 사용
        try:
            hd = primer3.calc_homodimer(seq, temp_c=oligo.tm)
            hd_dg = hd.dg / 1000.0 if hd.structure_found else 0.0
        except:
            hd_dg = 0.0
        
        # 3. 판정 (dG가 기준값보다 커야 통과. 예: -3.0 > -5.0)
        # criteria는 이제 Root QCCriteria이므로 hairpin_min_dg 접근 가능
        hp_pass = (hp_dg >= criteria.hairpin_min_dg)
        hd_pass = (hd_dg >= criteria.homodimer_min_dg)
        
        if not hp_pass: reasons.append("Hairpin")
        if not hd_pass: reasons.append("Homodimer")

        return {
            f"{prefix}_hairpin_dg": round(hp_dg, 2),
            f"{prefix}_homodimer_dg": round(hd_dg, 2),
            f"{prefix}_pass": hp_pass and hd_pass,
            f"{prefix}_fail_reason": "/".join(reasons)
        }

    def _check_hetero(self, s1: str, s2: str, criteria: QCCriteria, label: str) -> Dict[str, Any]:
        """두 서열 간의 Heterodimer 검사"""
        # ✅ 수정: snake_case 함수명 사용
        try:
            het = primer3.calc_heterodimer(s1, s2)
            het_dg = het.dg / 1000.0 if het.structure_found else 0.0
        except:
            het_dg = 0.0
        
        is_pass = (het_dg >= criteria.heterodimer_min_dg)
        
        return {
            f"{label}_dg": round(het_dg, 2),
            f"{label}_pass": is_pass
        }