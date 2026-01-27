from pydantic import BaseModel, Field
from typing import Optional

class ReferenceDetail(BaseModel):
    # YAML 예시: fasta: /path/to/hg19.fa
    fasta: str = Field(..., description="Path to Genome FASTA file")
    blast: str = Field(..., description="Path to BLAST database directory or prefix")

class ReferenceConfig(BaseModel):
    hg19: Optional[ReferenceDetail] = None
    hg38: Optional[ReferenceDetail] = None