import os

def create_file(path, content):
	print(os.path.dirname(path))
	os.makedirs(os.path.dirname(path), exist_ok=True)
	with open(path, "w") as f:
		f.write(content.strip() + "\n")
	print(f"Created: {path}")

# -------------------------------------------------------------------------
# 1. Packaging Configuration (Modern Standard)
# -------------------------------------------------------------------------
pyproject_toml = """
[build-system]
requires = ["setuptools>=61.0"]
build-backend = "setuptools.build_meta"

[project]
name = "pcr-pipeline"
version = "0.1.0"
description = "Bioinformatics PCR Primer Design Pipeline"
authors = [{name = "Bioinfo Developer", email = "dev@bioinfo.com"}]
requires-python = ">=3.9"
dependencies = [
	"pydantic>=2.0",
	"primer3-py>=2.0",
	"pyyaml"
]

[tool.setuptools.packages.find]
where = ["src"]
"""

# -------------------------------------------------------------------------
# 2. Schema (Pydantic Model) - Primer3 Option Mapping
# -------------------------------------------------------------------------
schema_code = """
from pydantic import BaseModel, Field, ConfigDict

class BasePrimer3Config(BaseModel):
	\"\"\"
	Primer3 옵션명(대문자)과 1:1 매핑되는 Config.
	YAML Preset에서 이 필드명 그대로 넘어오면 Designer가 자동으로 처리함.
	\"\"\"
	model_config = ConfigDict(extra='allow') # 정의되지 않은 추가 옵션 허용

	# 1. Task & General
	PRIMER_TASK: str = "generic"
	PRIMER_NUM_RETURN: int = 5
	PRIMER_PICK_INTERNAL_OLIGO: int = 0
		
	# 2. Tm Settings
	PRIMER_OPT_TM: float = 60.0
	PRIMER_MIN_TM: float = 55.0
	PRIMER_MAX_TM: float = 65.0
	PRIMER_PAIR_MAX_DIFF_TM: float = 5.0
		
	# 3. GC Settings
	PRIMER_OPT_GC_PERCENT: float = 50.0
	PRIMER_MIN_GC: float = 30.0
	PRIMER_MAX_GC: float = 70.0

	# 4. Helper Fields (Not direct Primer3 args, used for logic)
	amplicon_size_min: int = Field(80, description="Min Product Size")
	amplicon_size_max: int = Field(150, description="Max Product Size")
"""

# -------------------------------------------------------------------------
# 3. Loader (YAML Preset Loader)
# -------------------------------------------------------------------------
loader_code = """
import yaml
from pathlib import Path
from .schema import BasePrimer3Config

# 패키지 내부의 presets 폴더 경로 자동 추론
PRESET_DIR = Path(__file__).parent / "presets"

def load_preset(name: str = "default") -> BasePrimer3Config:
	yaml_path = PRESET_DIR / f"{name}.yaml"
		
	if not yaml_path.exists():
		# 파일이 없으면 기본값으로 생성 (혹은 에러 처리)
		if name == "default":
			return BasePrimer3Config()
		raise FileNotFoundError(f"Preset '{name}' not found at {yaml_path}")
		
	with open(yaml_path, "r") as f:
		data = yaml.safe_load(f) or {}
		
	# Pydantic을 통해 검증 및 객체 생성
	return BasePrimer3Config(**data)
"""

# -------------------------------------------------------------------------
# 4. YAML Presets
# -------------------------------------------------------------------------
yaml_default = """
# Default Settings
PRIMER_OPT_TM: 60.0
PRIMER_MIN_TM: 57.0
PRIMER_MAX_TM: 63.0
amplicon_size_min: 100
amplicon_size_max: 200
"""

yaml_high_gc = """
# Settings for High GC templates
PRIMER_OPT_TM: 64.0
PRIMER_MIN_TM: 60.0
PRIMER_MAX_TM: 68.0
PRIMER_OPT_GC_PERCENT: 60.0
PRIMER_MIN_GC: 50.0
PRIMER_MAX_GC: 80.0
PRIMER_SALT_DIVALENT: 2.5
"""

# -------------------------------------------------------------------------
# 5. Base Designer (Core Logic with Auto-Mapping)
# -------------------------------------------------------------------------
base_designer_code = """
from typing import Any, Dict, List, Optional
import primer3
from ..components import Amplicon, Primer
from ..config.schema import BasePrimer3Config

class BasePrimerDesigner:
	def __init__(
		self, 
		template_sequence: str,
		target_start: int,
		target_len: int,
		config: BasePrimer3Config,
		overrides: Optional[Dict[str, Any]] = None
	):
		self.seq = template_sequence
		self.target = [target_start, target_len]
		
		# 1. Config 적용 (Overrides가 있으면 덮어쓰기)
		if overrides:
			self.cfg = config.model_copy(update=overrides)
		else:
			self.cfg = config
			
		self.global_args = {}
		self.seq_args = {}
		self._init_args()

	def _init_args(self):
		# A. Sequence Args
		self.seq_args = {
			"SEQUENCE_ID": "demo_run",
			"SEQUENCE_TEMPLATE": self.seq,
			"SEQUENCE_TARGET": self.target
		}

		# B. Global Args (Auto-Mapping)
		# Config 객체를 dict로 변환 후, 'PRIMER_'로 시작하는 키만 추출
		for key, value in self.cfg.model_dump().items():
			if key.startswith("PRIMER_"):
				self.global_args[key] = value
				
		# C. 로직 처리 필요한 필드 (Range 등)
		self.global_args["PRIMER_PRODUCT_SIZE_RANGE"] = [
			[self.cfg.amplicon_size_min, self.cfg.amplicon_size_max]
		]

	def design(self) -> List[Amplicon]:
		# Primer3 실행
		res = primer3.bindings.designPrimers(self.seq_args, self.global_args)
		
		# 결과 파싱 (간소화된 버전)
		amplicons = []
		count = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))
		for i in range(count):
			fwd_seq = res.get(f"PRIMER_LEFT_{i}_SEQUENCE")
			rev_seq = res.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
			if fwd_seq and rev_seq:
				# 여기서 Components의 객체를 생성
				amplicons.append(Amplicon(
					fwd_seq=fwd_seq, 
					rev_seq=rev_seq,
					tm=res.get(f"PRIMER_LEFT_{i}_TM")
				))
		return amplicons
"""

# -------------------------------------------------------------------------
# 6. Components (Simple Data Classes)
# -------------------------------------------------------------------------
components_code = """
from dataclasses import dataclass

@dataclass
class Primer:
	sequence: str

@dataclass
class Amplicon:
	fwd_seq: str
	rev_seq: str
	tm: float
"""

# -------------------------------------------------------------------------
# 7. Main Script (Local Execution Test)
# -------------------------------------------------------------------------
run_local_code = """
from pcr.config.loader import load_preset
from pcr.designers.base import BasePrimerDesigner

# 1. Preset 로드 (YAML 파일 이름만 지정)
# config/presets/high_gc.yaml 을 읽어옵니다.
config = load_preset("high_gc")

print(f"Loaded Config: Target Tm = {config.PRIMER_OPT_TM}")

# 2. 실행 (필요하다면 코드 레벨에서 추가 오버라이드 가능)
seq = "ATGC" * 100 # Dummy Sequence
designer = BasePrimerDesigner(
	template_sequence=seq, 
	target_start=50, 
	target_len=20, 
	config=config,
	overrides={"PRIMER_NUM_RETURN": 1} # 코드에서 즉석 변경 예시
)

# 3. 디자인
results = designer.design()
print(f"Designed {len(results)} amplicons.")
for amp in results:
	print(f" - Fwd: {amp.fwd_seq} / Tm: {amp.tm}")
"""

# -------------------------------------------------------------------------
# Generate Files
# -------------------------------------------------------------------------
# Package Root: src/pcr
base_dir = "src/pcr"

files = {
	"/storage/home/jhkim/scripts/Task/pcr_designer/pcrr/pyproject.toml": pyproject_toml,
	f"{base_dir}/__init__.py": "",
	f"{base_dir}/components/__init__.py": components_code,  # Simplified for demo
	f"{base_dir}/config/__init__.py": "",
	f"{base_dir}/config/schema.py": schema_code,
	f"{base_dir}/config/loader.py": loader_code,
	f"{base_dir}/config/presets/default.yaml": yaml_default,
	f"{base_dir}/config/presets/high_gc.yaml": yaml_high_gc,
	f"{base_dir}/designers/__init__.py": "",
	f"{base_dir}/designers/base.py": base_designer_code,
	"/storage/home/jhkim/scripts/Task/pcr_designer/pcrr/run_local.py": run_local_code
}

for path, content in files.items():
	create_file(path, content)

print("\n✅ Package structure created successfully!")
print("👉 Next Step: Run 'pip install -e .' to install config and core logic as a package.")