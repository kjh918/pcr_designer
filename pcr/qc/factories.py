# pcr/qc/factories.py
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from pcr.components import Primer, Amplicon


def _safe_setattr(obj: Any, key: str, value: Any) -> None:
    try:
        setattr(obj, key, value)
    except Exception:
        pass


def _mk_primer_with_binding(
    *,
    template: str,
    seq: str,
    strand: str,
    primer_type: str,
    target_start_index: int,
    target_end_index: int,
    binding_start_index: Optional[int] = None,
    binding_end_index: Optional[int] = None,
) -> Primer:
    """
    ✅ CHANGED: binding_start_index / binding_end_index를 Primer에 주입
    """
    kwargs: Dict[str, Any] = dict(
        template_sequence=template,
        reference_template_sequence=template,
        sequence=seq,
        strand=strand,
        primer_type=primer_type,
        target_start_index=target_start_index,
        target_end_index=target_end_index,
    )

    # ✅ ADDED
    if binding_start_index is not None:
        kwargs["binding_start_index"] = int(binding_start_index)
    if binding_end_index is not None:
        kwargs["binding_end_index"] = int(binding_end_index)

    try:
        p = Primer(**kwargs)  # type: ignore[arg-type]
    except TypeError:
        # ✅ ADDED: Primer 생성자가 해당 필드를 못 받으면 setattr로 부착
        base_kwargs = dict(
            template_sequence=template,
            reference_template_sequence=template,
            sequence=seq,
            strand=strand,
            primer_type=primer_type,
            target_start_index=target_start_index,
            target_end_index=target_end_index,
        )
        p = Primer(**base_kwargs)  # type: ignore[arg-type]
        _safe_setattr(p, "binding_start_index", binding_start_index)
        _safe_setattr(p, "binding_end_index", binding_end_index)

    return p


# ============================
# ✅ NEW: filtered_amplicons(여러개) -> Amplicon list
# ============================
def make_amplicons_for_qc_from_blast(
    *,
    blast_result: Dict[str, Any],
    forward_seq: str,
    reverse_seq: str,
    probe_seq: str | None = None,
) -> List[Amplicon]:
    """
    ✅ CHANGED: blast_qc_for_primer_pair() 결과 전체를 받아서
    result.filtered_amplicons만 Amplicon 리스트로 만들어 반환
    """
    if not blast_result or blast_result.get("blast_error") is True:
        return []

    result = blast_result.get("result") or {}
    filtered_amplicons = result.get("filtered_amplicons") or []
    if not isinstance(filtered_amplicons, list):
        return []

    amps: List[Amplicon] = []
    for amp_dict in filtered_amplicons:
        try:
            amps.append(
                make_amplicon_for_qc_from_blast_amplicon(
                    blast_amplicon=amp_dict,
                    forward_seq=forward_seq,
                    reverse_seq=reverse_seq,
                    probe_seq=probe_seq,
                )
            )
        except Exception:
            # 하나 실패해도 나머지는 진행
            continue

    return amps


# ============================
# ✅ NEW: filtered_amplicons의 "원소 1개" -> Amplicon 1개
# ============================
def make_amplicon_for_qc_from_blast_amplicon(
    *,
    blast_amplicon: Dict[str, Any],
    forward_seq: str,
    reverse_seq: str,
    probe_seq: str | None = None,
) -> Amplicon:
    """
    blast_amplicon 구조(너 blast.py에서 만든 one dict):
      {
        "PASS": bool,
        "reference": {"start_1b","end_1b","seq","forward_full_bind_1b":{start,end},"reverse_full_bind_1b":{start,end}},
        "amplicon": {"start_1b","end_1b",...},
        ...
      }
    """
    forward_seq = forward_seq.replace(" ", "").strip().upper()
    reverse_seq = reverse_seq.replace(" ", "").strip().upper()
    probe_seq = probe_seq.replace(" ", "").strip().upper() if probe_seq else None

    ref = blast_amplicon.get("reference") or {}
    amp = blast_amplicon.get("amplicon") or {}

    template = (ref.get("seq") or "").strip().upper()
    if not template:
        # ✅ FIX: 예전 에러가 여기서 발생하던 것
        raise ValueError("No reference template seq in blast_amplicon.reference.seq")

    ref_start_1b = int(ref.get("start_1b"))
    # ref_end_1b = int(ref.get("end_1b"))  # 필요하면 사용

    # ✅ target = primer 사이 inner amplicon 영역(0-based, inclusive)
    amp_start_1b = int(amp.get("start_1b"))
    amp_end_1b = int(amp.get("end_1b"))

    if amp_start_1b <= amp_end_1b:
        target_start_index = max(0, amp_start_1b - ref_start_1b)
        target_end_index = min(len(template) - 1, amp_end_1b - ref_start_1b)
    else:
        # 사이 영역이 없으면 전체를 target으로
        target_start_index = 0
        target_end_index = len(template) - 1

    # ✅ binding index (0-based, inclusive) = full_bind_1b를 template 기준으로 변환
    f_full = ref.get("forward_full_bind_1b") or {}
    r_full = ref.get("reverse_full_bind_1b") or {}

    f_bind_s_1b = int(f_full.get("start"))
    f_bind_e_1b = int(f_full.get("end"))
    r_bind_s_1b = int(r_full.get("start"))
    r_bind_e_1b = int(r_full.get("end"))

    f_bind_start_index = max(0, f_bind_s_1b - ref_start_1b)
    f_bind_end_index = min(len(template) - 1, f_bind_e_1b - ref_start_1b)
    r_bind_start_index = max(0, r_bind_s_1b - ref_start_1b)
    r_bind_end_index = min(len(template) - 1, r_bind_e_1b - ref_start_1b)

    f_primer = _mk_primer_with_binding(
        template=template,
        seq=forward_seq,
        strand="forward",
        primer_type="forward",
        target_start_index=target_start_index,
        target_end_index=target_end_index,
        binding_start_index=f_bind_start_index,
        binding_end_index=f_bind_end_index,
    )

    r_primer = _mk_primer_with_binding(
        template=template,
        seq=reverse_seq,
        strand="reverse",
        primer_type="reverse",
        target_start_index=target_start_index,
        target_end_index=target_end_index,
        binding_start_index=r_bind_start_index,
        binding_end_index=r_bind_end_index,
    )

    probe_primer: Primer | None = None
    if probe_seq:
        probe_primer = _mk_primer_with_binding(
            template=template,
            seq=probe_seq,
            strand="reverse",
            primer_type="probe",
            target_start_index=target_start_index,
            target_end_index=target_end_index,
            binding_start_index=None,
            binding_end_index=None,
        )

    return Amplicon(
        template_sequence=template,
        reference_template_sequence=template,
        target_start_index=target_start_index,
        target_end_index=target_end_index,
        forward_primer=f_primer,
        reverse_primer=r_primer,
        probe=probe_primer,
    )
