from pathlib import Path
from typing import Dict
import yaml
from pydantic import BaseModel

BASE_DIR = Path(__file__).resolve().parents[1]   # project root (app, primer 있는 곳)
CONFIG_FILE = BASE_DIR / "config" / "references.yaml"


class ReferenceConfig(BaseModel):
    name: str
    fasta: Path


class Settings(BaseModel):
    references: Dict[str, ReferenceConfig]


def load_settings() -> Settings:
    with open(CONFIG_FILE, "r") as f:
        raw = yaml.safe_load(f)

    refs: Dict[str, ReferenceConfig] = {}
    for name, cfg in raw.get("references", {}).items():
        refs[name] = ReferenceConfig(
            name=name,
            fasta=Path(cfg["fasta"]),
        )

    return Settings(references=refs)


settings = load_settings()
