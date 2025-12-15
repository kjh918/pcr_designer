from pathlib import Path
from pydantic import BaseModel

class ReferenceConfig(BaseModel):
    fasta: Path
    blast: Path
