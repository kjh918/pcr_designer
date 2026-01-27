from pydantic import BaseModel
from .references import ReferenceConfig
from .pcr import PCRParams
from .qc import QCParams

class AppConfig(BaseModel):
    references: ReferenceConfig
    pcr_params: PCRParams
    qc_params: QCParams