from typing import Dict, Any
import primer3
from ..config.schema.qc import QCParams, BaseQCCriteria
from ..components import Amplicon, Primer
from .types import QCResult

class ThermoChecker:
    def __init__(self, qc_params: QCParams):
        self.primer_criteria = qc_params.get_primer_criteria()
        self.probe_criteria = qc_params.get_probe_criteria()

    def check(self, amplicon: Amplicon) -> QCResult:
        """물리적 성질 계산 후 QCResult 반환 (객체 수정 X)"""
        details = {}
        all_pass = True

        # 1. Forward
        f_res = self._check_oligo(amplicon.forward_primer, self.primer_criteria, "fwd")
        details.update(f_res)
        if not f_res["fwd_pass"]: all_pass = False

        # 2. Reverse
        r_res = self._check_oligo(amplicon.reverse_primer, self.primer_criteria, "rev")
        details.update(r_res)
        if not r_res["rev_pass"]: all_pass = False

        # 3. Heterodimer
        fr_res = self._check_hetero(amplicon.forward_primer.sequence, 
                                    amplicon.reverse_primer.sequence, 
                                    self.primer_criteria, "hetero_fr")
        details.update(fr_res)
        if not fr_res["hetero_fr_pass"]: all_pass = False

        # 4. Probe (Optional)
        if amplicon.probe:
            p_res = self._check_oligo(amplicon.probe, self.probe_criteria, "probe")
            details.update(p_res)
            if not p_res["probe_pass"]: all_pass = False
            # ... (Probe Hetero 생략, 필요시 추가) ...

        return QCResult(passed=all_pass, data=details)

    def _check_oligo(self, oligo: Primer, criteria: BaseQCCriteria, prefix: str) -> Dict[str, Any]:
        # (기존 로직 유지하되 리턴값만 dict)
        hp = primer3.calc_hairpin(oligo.sequence)
        hp_dg = hp.dg / 1000.0 if hp.structure_found else 0.0
        
        hd = primer3.calc_homodimer(oligo.sequence)
        hd_dg = hd.dg / 1000.0 if hd.structure_found else 0.0
        
        pass_flag = (hp_dg >= criteria.hairpin_min_dg) and (hd_dg >= criteria.homodimer_min_dg)
        
        return {
            f"{prefix}_hairpin_dg": hp_dg,
            f"{prefix}_homodimer_dg": hd_dg,
            f"{prefix}_pass": pass_flag
        }

    def _check_hetero(self, s1: str, s2: str, criteria: BaseQCCriteria, label: str) -> Dict[str, Any]:
        het = primer3.calc_heterodimer(s1, s2)
        het_dg = het.dg / 1000.0 if het.structure_found else 0.0
        return {
            f"{label}_dg": het_dg,
            f"{label}_pass": het_dg >= criteria.heterodimer_min_dg
        }