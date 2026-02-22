#!/bin/bash
# init_project.sh
# 실행: bash init_project.sh [프로젝트 루트 경로]
# 예시: bash init_project.sh /home/user/myproject

ROOT=${1:-.}

echo "📁 Creating project structure at: $ROOT"

# ─────────────────────────────────────────
# 디렉토리 생성
# ─────────────────────────────────────────
mkdir -p \
  $ROOT/pcr/components \
  $ROOT/pcr/designers/base \
  $ROOT/pcr/designers/qpcr \
  $ROOT/pcr/designers/as_pcr \
  $ROOT/pcr/designers/ms_pcr \
  $ROOT/pcr/qc \
  $ROOT/pcr/utils \
  $ROOT/pcr/config/schema \
  $ROOT/scripts \
  $ROOT/test

# ─────────────────────────────────────────
# __init__.py
# ─────────────────────────────────────────
touch $ROOT/pcr/__init__.py
touch $ROOT/pcr/components/__init__.py
touch $ROOT/pcr/designers/__init__.py
touch $ROOT/pcr/designers/base/__init__.py
touch $ROOT/pcr/designers/qpcr/__init__.py
touch $ROOT/pcr/designers/as_pcr/__init__.py
touch $ROOT/pcr/designers/ms_pcr/__init__.py
touch $ROOT/pcr/qc/__init__.py
touch $ROOT/pcr/utils/__init__.py
touch $ROOT/pcr/config/__init__.py
touch $ROOT/pcr/config/schema/__init__.py

# ─────────────────────────────────────────
# components/
# ─────────────────────────────────────────
touch $ROOT/pcr/components/region.py
touch $ROOT/pcr/components/primer.py
touch $ROOT/pcr/components/amplicon.py

# ─────────────────────────────────────────
# designers/base/
# ─────────────────────────────────────────
touch $ROOT/pcr/designers/base/designer.py
touch $ROOT/pcr/designers/base/schema.py

# ─────────────────────────────────────────
# designers/qpcr/
# ─────────────────────────────────────────
touch $ROOT/pcr/designers/qpcr/designer.py
touch $ROOT/pcr/designers/qpcr/schema.py

# ─────────────────────────────────────────
# designers/as_pcr/
# ─────────────────────────────────────────
touch $ROOT/pcr/designers/as_pcr/designer.py
touch $ROOT/pcr/designers/as_pcr/schema.py

# ─────────────────────────────────────────
# designers/ms_pcr/
# ─────────────────────────────────────────
touch $ROOT/pcr/designers/ms_pcr/designer.py
touch $ROOT/pcr/designers/ms_pcr/schema.py

# ─────────────────────────────────────────
# qc/
# ─────────────────────────────────────────
touch $ROOT/pcr/qc/executor.py
touch $ROOT/pcr/qc/blast.py
touch $ROOT/pcr/qc/thermo.py
touch $ROOT/pcr/qc/ispcr.py

# ─────────────────────────────────────────
# utils/
# ─────────────────────────────────────────
touch $ROOT/pcr/utils/ranker.py

# ─────────────────────────────────────────
# config/
# ─────────────────────────────────────────
touch $ROOT/pcr/config/loader.py
touch $ROOT/pcr/config/system.yaml
touch $ROOT/pcr/config/base_pcr.yaml
touch $ROOT/pcr/config/schema/app.py
touch $ROOT/pcr/config/schema/pcr.py
touch $ROOT/pcr/config/schema/qc.py
touch $ROOT/pcr/config/schema/references.py

# ─────────────────────────────────────────
# factory / scripts
# ─────────────────────────────────────────
touch $ROOT/pcr/factory.py
touch $ROOT/scripts/run_pipeline.py
touch $ROOT/scripts/test.py

# ─────────────────────────────────────────
# 결과 출력
# ─────────────────────────────────────────
echo ""
echo "✅ Done. Structure:"
find $ROOT -not -path "*/__pycache__/*" | sort | \
  awk '{
    n = split($0, a, "/")
    indent = ""
    for (i=2; i<n; i++) indent = indent "│   "
    if (n > 1) {
      if ($0 ~ /\/$/ || system("test -d "$0) == 0)
        print indent "├── " a[n] "/"
      else
        print indent "├── " a[n]
    } else {
      print a[n]
    }
  }'