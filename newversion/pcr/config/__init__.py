#from .loader import ConfigLoader
from .schema.root import AppConfig, BaseDesignInput, BaseDesignOutput
from .schema.pcr import PCRParams, PrimerKwargs, ProbeKwargs
# QCParams 대신 실제 클래스 이름인 QCCriteria, QCToolsConfig를 가져옵니다.
from .schema.qc import QCCriteria, QCToolsConfig, PrimerQCCriteria, ProbeQCCriteria

__all__ = [
    "AppConfig",
    "BaseDesignInput",
    "BaseDesignOutput",
    "PCRParams",
    "QCCriteria",      # 수정
    "QCToolsConfig",   # 수정
    "PrimerKwargs",
    "ProbeKwargs",
    "PrimerQCCriteria",
    "ProbeQCCriteria"
]