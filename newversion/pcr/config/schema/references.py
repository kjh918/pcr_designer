from pydantic import BaseModel, Field

class ReferenceConfig(BaseModel):
    """
    system.yaml의 references 섹션 매핑
    """
    # 수정: YAML 키 이름(fasta_path)과 일치시킴
    fasta_path: str = Field(..., alias="fasta_path") 
    
    # 수정: YAML 키 이름(blast_db_path)과 일치시킴
    blast_db_path: str = Field(..., alias="blast_db_path")

    class Config:
        populate_by_name = True