import primer3
from typing import Dict, Any
from ..config.schema.qc import QCCriteria
from ..components.amplicon import Amplicon, Primer

class ThermoChecker:
    """
    물리적 성질(Hairpin, Dimer dG)을 검사하는 클래스.
    Primer3-py를 사용하여 dG(Gibbs Free Energy, kcal/mol)를 산출합니다.
    """
    def __init__(self, criteria: QCCriteria):
        self.criteria = criteria

    def check(self, amplicon: Amplicon) -> Dict[str, Any]:
        """
        Amplicon 내의 모든 올리고(F, R, P)에 대해 물리적 안정성을 검사합니다.
        """
        details = {}
        all_pass = True
        fail_reasons = []

        # ---------------------------------------------------------
        # 1. Individual Oligo Check (Hairpin & Homodimer)
        # ---------------------------------------------------------
        
        # (1) Forward Primer
        f_res = self._check_oligo(amplicon.forward, "fwd")
        details.update(f_res)
        if not f_res["fwd_pass"]:
            all_pass = False
            fail_reasons.append(f"Fwd_Thermo({f_res.get('fwd_fail_reason')})")

        # (2) Reverse Primer
        r_res = self._check_oligo(amplicon.reverse, "rev")
        details.update(r_res)
        if not r_res["rev_pass"]:
            all_pass = False
            fail_reasons.append(f"Rev_Thermo({r_res.get('rev_fail_reason')})")

        # (3) Probe (Optional)
        if amplicon.probe:
            p_res = self._check_oligo(amplicon.probe, "probe")
            details.update(p_res)
            if not p_res["probe_pass"]:
                all_pass = False
                fail_reasons.append(f"Probe_Thermo({p_res.get('probe_fail_reason')})")
        else:
            # Probe가 없으면 기본값 채움 (리포트 에러 방지)
            details.update({
                "probe_hairpin_dg": 0.0,
                "probe_homodimer_dg": 0.0,
                "probe_pass": True
            })

        # ---------------------------------------------------------
        # 2. Heterodimer Check (Interactions)
        # ---------------------------------------------------------
        
        # (1) Forward + Reverse (Primer Dimer) - 필수
        fr_res = self._check_hetero(
            amplicon.forward.sequence, 
            amplicon.reverse.sequence, 
            "hetero_fr"
        )
        details.update(fr_res)
        if not fr_res["hetero_fr_pass"]:
            all_pass = False
            fail_reasons.append("Hetero_FR")

        # (2) Probe Interactions (FP, RP) - 선택
        if amplicon.probe:
            # Forward + Probe
            fp_res = self._check_hetero(
                amplicon.forward.sequence, 
                amplicon.probe.sequence, 
                "hetero_fp"
            )
            # Reverse + Probe
            rp_res = self._check_hetero(
                amplicon.reverse.sequence, 
                amplicon.probe.sequence, 
                "hetero_rp"
            )
            
            details.update(fp_res)
            details.update(rp_res)

            if not (fp_res["hetero_fp_pass"] and rp_res["hetero_rp_pass"]):
                all_pass = False
                fail_reasons.append("Hetero_Probe")
        else:
            # Probe 없으면 0.0 처리
            details.update({
                "hetero_fp_dg": 0.0, "hetero_fp_pass": True,
                "hetero_rp_dg": 0.0, "hetero_rp_pass": True
            })

        # ---------------------------------------------------------
        # 3. Final Return
        # ---------------------------------------------------------
        return {
            "passed": all_pass,
            "fail_reason": ", ".join(fail_reasons) if fail_reasons else "",
            "data": details  # 이 딕셔너리가 summarize로 전달됨
        }

    def _check_oligo(self, oligo: Primer, prefix: str) -> Dict[str, Any]:
        """개별 올리고뉴클레오티드의 Hairpin 및 Homodimer 검사"""
        seq = oligo.sequence
        tm = oligo.tm
        reasons = []
        
        # 기준값 (kcal/mol, 예: -5.0)
        # dG 값이 이보다 커야 통과 (예: -3.0 > -5.0)
        min_hp = self.criteria.hairpin_min_dg
        min_hd = self.criteria.homodimer_min_dg

        # 1. Hairpin
        try:
            hp = primer3.calc_hairpin(seq, temp_c=tm)
            # primer3는 cal/mol 반환 -> kcal/mol로 변환 (/1000)
            hp_dg = hp.dg / 1000.0 if hp.structure_found else 0.0
        except:
            hp_dg = 0.0

        # 2. Homodimer
        try:
            hd = primer3.calc_homodimer(seq, temp_c=tm)
            hd_dg = hd.dg / 1000.0 if hd.structure_found else 0.0
        except:
            hd_dg = 0.0

        # 3. 판정
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
        """두 서열 간의 Heterodimer 검사"""
        min_het = self.criteria.heterodimer_min_dg
        
        try:
            het = primer3.calc_heterodimer(s1, s2)
            het_dg = het.dg / 1000.0 if het.structure_found else 0.0
        except:
            het_dg = 0.0
        
        is_pass = (het_dg >= min_het)
        
        return {
            f"{label}_dg": round(het_dg, 2),
            f"{label}_pass": is_pass
        }