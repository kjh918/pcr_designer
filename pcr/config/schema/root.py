from typing import Dict
from pydantic import BaseModel, Field, model_validator

from config.schema.references import ReferenceConfig
from config.schema.pcr import PCRParams
from config.schema.qc import QCParams

class Settings(BaseModel):
    references: Dict[str, ReferenceConfig] = Field(default_factory=dict)
    pcr_params: PCRParams = Field(default_factory=PCRParams)
    qc_params: QCParams

    @model_validator(mode="after")
    def _validate_refs(self) -> "Settings":
        if not self.references:
            raise ValueError("No references configured (references: ...).")
        return self
