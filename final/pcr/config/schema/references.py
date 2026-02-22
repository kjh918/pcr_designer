from pydantic import BaseModel
from typing import Optional

class SystemConfig(BaseModel):
    blastn_path: str = "blastn"
    makeblastdb_path: Optional[str] = None

class ReferenceConfig(BaseModel):
    fasta_path: str
    blast_db_path: str