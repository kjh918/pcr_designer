# primer/_get_cfg.py
from typing import Dict, Optional, Tuple
from pyfaidx import Fasta

try:
    from config.settings import settings
except ImportError:
    settings = None  # 테스트용 등 settings 없는 환경 대비

# ref_name(hg19/hg38) 별 Fasta 핸들 캐시
_FASTA_CACHE: Dict[str, Fasta] = {}


def get_fasta_handle(ref_name: str) -> Fasta:
    """
    config.settings.references[ref_name].fasta 를 이용해서
    pyfaidx.Fasta 객체를 가져온다 (캐시 사용).
    """
    if settings is None:
        raise RuntimeError("config.settings 를 불러올 수 없습니다. settings 설정을 확인하세요.")

    if ref_name not in settings.references:
        raise ValueError(f"Unknown reference name: {ref_name}")

    if ref_name not in _FASTA_CACHE:
        ref_cfg = settings.references[ref_name]
        _FASTA_CACHE[ref_name] = str(ref_cfg.fasta)

    return _FASTA_CACHE[ref_name]


def get_pcr_params_with_override(
    min_amplicon_length: Optional[int],
    max_amplicon_length: Optional[int],
    n_probes: Optional[int],
    n_primers: Optional[int],
    bisulfite: Optional[bool],
) -> Tuple[int, int, int, int, bool]:
    """
    ⚙️ 웹에서 넘어온 값이 있으면 그 값을 우선으로 사용하고,
    None 인 항목만 config.settings.pcr_params 값으로 채운다.

    반환: (min_amplicon_length, max_amplicon_length, n_probes, n_primers, bisulfite)
    """
    if settings is None:
        # settings 없이 돌리는 환경이면, 모든 값을 직접 넘기도록 강제
        missing = [
            name
            for name, val in [
                ("min_amplicon_length", min_amplicon_length),
                ("max_amplicon_length", max_amplicon_length),
                ("n_probes", n_probes),
                ("n_primers", n_primers),
                ("bisulfite", bisulfite),
            ]
            if val is None
        ]
        if missing:
            raise RuntimeError(
                f"PCR 파라미터 {missing} 가 None 이고, config.settings 를 사용할 수 없습니다. "
                "테스트/CLI 환경이라면 인자를 모두 직접 넘겨주세요."
            )
        return (
            min_amplicon_length,  # type: ignore
            max_amplicon_length,  # type: ignore
            n_probes,             # type: ignore
            n_primers,            # type: ignore
            bisulfite,           # type: ignore
        )

    pcr_cfg = settings.pcr_params

    # ▶ 웹에서 온 값이 None 이면 config 값으로 채우기
    if min_amplicon_length is None:
        min_amplicon_length = pcr_cfg.primer_kwargs.min_amplicon_length
    if max_amplicon_length is None:
        max_amplicon_length = pcr_cfg.primer_kwargs.max_amplicon_length
    if n_primers is None:
        n_primers = pcr_cfg.primer_kwargs.n_primers
    if n_probes is None:
        n_probes = pcr_cfg.probe_kwargs.n_probes
    if bisulfite is None:
        bisulfite = pcr_cfg.bisulfite.run

    return (
        min_amplicon_length,
        max_amplicon_length,
        n_probes,
        n_primers,
        bisulfite,
    )
