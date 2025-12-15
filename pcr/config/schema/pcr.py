from typing import Any, Dict
from pydantic import BaseModel, Field

class PrimerKwargs(BaseModel):
    min_amplicon_length: int = 80
    max_amplicon_length: int = 120
    n_primers: int = 100
    opt_length: int = 25
    min_length: int = 20
    max_length: int = 30
    opt_tm: float = 60
    min_tm: float = 55
    max_tm: float = 65
    opt_gc: float = 50
    min_gc: float = 35
    max_gc: float = 65
    primer3_global_args: Dict[str, Any] = Field(default_factory=dict)

class ProbeKwargs(BaseModel):
    n_probes: int = 100
    min_primer_probe_tm_diff: float = 6
    max_primer_probe_tm_diff: float = 8
    opt_length: int = 25
    min_length: int = 20
    max_length: int = 30
    opt_tm: float = 60
    min_tm: float = 60
    max_tm: float = 65
    opt_gc: float = 50
    min_gc: float = 35
    max_gc: float = 65
    primer3_global_args: Dict[str, Any] = Field(default_factory=dict)

class BisulfiteConfig(BaseModel):
    run: bool = False
    cpg_default: str = "methyl"

class PCRParams(BaseModel):
    primer_kwargs: PrimerKwargs = Field(default_factory=PrimerKwargs)
    probe_kwargs: ProbeKwargs = Field(default_factory=ProbeKwargs)
    bisulfite: BisulfiteConfig = Field(default_factory=BisulfiteConfig)
