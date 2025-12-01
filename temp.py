"""
===========================================================
  qPCR Full Pipeline (Single Script Consolidation)
  Author: ChatGPT
===========================================================
  
  이 스크립트 하나로 다음 기능을 모두 포함한다:

  1) PrimerDesigner  
     - primer-only 모드  
     - probe 모드  
     - primer3 기반 Tm/GC/self-dimer/hairpin/hetero-dimer 계산  
     - primer3 글로벌 옵션 자동 설정  

  2) BLAST specificity check
     - blastn-short 사용  
     - probe/primer off-target 정량  

  3) Amplicon scoring
     - primer3 thermo penalties  
     - BLAST off-target penalties  
     - 최종 score 기반 ranking  

  4) FastAPI
     - /api/qpcr → qPCR 분석 실행  
     - 결과 파일(JSON) 자동 저장  
     - /download/json/{job_id}, /download/tsv/{job_id} 제공  

  5) 디렉토리 자동 생성
     - data/results/ 아래 job_id.json 저장  
===========================================================
"""

############################################################
# 0. Imports
############################################################

import os
import json
import uuid
import tempfile
import subprocess
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Optional

from fastapi import FastAPI, HTTPException
from fastapi.responses import JSONResponse, FileResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

import primer3

############################################################
# 1. Primer & Amplicon 데이터 구조
############################################################

class Primer:
    """단일 primer 객체"""
    def __init__(self, sequence: str, strand: str, primer_type: str):
        self.sequence = sequence
        self.strand = strand
        self.primer_type = primer_type


class Amplicon:
    """한 primer set으로 만들어진 amplicon 객체"""
    def __init__(self, forward: Primer, reverse: Primer, probe: Optional[Primer]):
        self.forward_primer = forward
        self.reverse_primer = reverse
        self.probe = probe

        self.final_score = None   # scoring 단계에서 평가됨
        self.f_self_dimer = None
        self.r_self_dimer = None
        self.hetero_dimer = None
        self.f_offtarget = None
        self.r_offtarget = None
        self.probe_offtarget = None


############################################################
# 2. Primer3 Thermo 체크 (self-dimer, hetero-dimer, hairpin)
############################################################

def thermo_check(seq_f: str, seq_r: Optional[str] = None):
    """
    seq_f, seq_r 에 대하여:
      - forward self-dimer
      - reverse self-dimer
      - hairpin
      - hetero-dimer (forward vs reverse)
    """
    result = {
        "f_self_dimer": primer3.calcHomodimer(seq_f).dg,
        "f_hairpin": primer3.calcHairpin(seq_f).dg
    }

    if seq_r:
        result["r_self_dimer"] = primer3.calcHomodimer(seq_r).dg
        result["r_hairpin"] = primer3.calcHairpin(seq_r).dg
        result["hetero_dimer"] = primer3.calcHeterodimer(seq_f, seq_r).dg

    return result


############################################################
# 3. BLAST specificity scoring 함수
############################################################

def blast_offtargets(seq: str, blast_db: str) -> int:
    """
    blastn-short 사용하여 off-target 개수 계수
    - 자기 자신 매칭 1개는 제거
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        qpath = Path(tmpdir) / "q.fa"
        qpath.write_text(f">q\n{seq}\n")

        cmd = [
            "blastn", "-task", "blastn-short",
            "-query", str(qpath),
            "-db", blast_db,
            "-outfmt", "6"
        ]

        p = subprocess.run(cmd, capture_output=True, text=True)
        hits = [line for line in p.stdout.split("\n") if line.strip()]

        return max(0, len(hits) - 1)  # 본인 포함 1개 빼기


############################################################
# 4. PrimerDesigner: primer-only, probe 모드 모두 지원
############################################################

class PrimerDesigner:
    """
    mode = 'primer' 또는 'probe'
    """
    def __init__(self, sequence: str, mode="primer"):
        self.sequence = sequence.upper()
        self.mode = mode
        self.amplicon_list = []

    def design(self):
        """
        primer3.bindings.designPrimers 를 실행하고 Amplicon 리스트를 생성
        """
        seq_args = {
            "SEQUENCE_ID": "TARGET",
            "SEQUENCE_TEMPLATE": self.sequence
        }
        global_args = {
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": int(self.mode == "probe"),
            "PRIMER_NUM_RETURN": 10,
            "PRIMER_OPT_SIZE": 22,
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_MAX_SIZE": 30,
            "PRIMER_PRODUCT_SIZE_RANGE": [[80, 180]],
        }

        res = primer3.bindings.designPrimers(seq_args, global_args)

        total = res.get("PRIMER_PAIR_NUM_RETURNED", 0)

        amps = []
        for i in range(total):
            left = res.get(f"PRIMER_LEFT_{i}_SEQUENCE")
            right = res.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
            internal = res.get(f"PRIMER_INTERNAL_{i}_SEQUENCE")

            f = Primer(left, "forward", "forward") if left else None
            r = Primer(right, "reverse", "reverse") if right else None
            p = Primer(internal, "forward", "probe") if internal else None

            amps.append(Amplicon(f, r, p))

        self.amplicon_list = amps
        return amps


############################################################
# 5. Amplicon Scoring (thermo + BLAST)
############################################################

def score_amplicon(amp: Amplicon, blast_db: str):
    """최종 score = thermo penalty + BLAST penalty"""
    f = amp.forward_primer.sequence
    r = amp.reverse_primer.sequence
    p = amp.probe.sequence if amp.probe else None

    # Thermo
    t = thermo_check(f, r)
    amp.f_self_dimer = t["f_self_dimer"]
    amp.r_self_dimer = t.get("r_self_dimer")
    amp.hetero_dimer = t.get("hetero_dimer")

    thermo_penalty = (
        abs(t["f_self_dimer"]) +
        abs(t["f_hairpin"]) +
        abs(t.get("r_self_dimer", 0)) +
        abs(t.get("r_hairpin", 0)) +
        abs(t.get("hetero_dimer", 0))
    )

    # BLAST off-target
    f_ot = blast_offtargets(f, blast_db)
    r_ot = blast_offtargets(r, blast_db)
    p_ot = blast_offtargets(p, blast_db) if p else 0

    amp.f_offtarget = f_ot
    amp.r_offtarget = r_ot
    amp.probe_offtarget = p_ot

    blast_penalty = 3*f_ot + 3*r_ot + 5*p_ot

    # 최종 점수
    score = thermo_penalty + blast_penalty
    amp.final_score = score
    return score


############################################################
# 6. FastAPI 서버 + 결과 저장 & 다운로드
############################################################

app = FastAPI()

# templates, static
app.mount("/static", StaticFiles(directory="app/static"), name="static")
templates = Jinja2Templates(directory="app/templates")

# 결과 폴더
BASE = Path(__file__).resolve().parent
RESULTS = BASE / "data" / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


############################################################
# Request Schema (간단 버전)
############################################################

from pydantic import BaseModel

class QpcrRequest(BaseModel):
    sequence: str
    mode: str = "probe"   # "primer" or "probe"
    blast_db: str         # e.g., "/path/to/hg38/blastdb_ref"


############################################################
# API: qPCR 분석 실행
############################################################

@app.post("/api/qpcr")
def run_qpcr(req: QpcrRequest):

    designer = PrimerDesigner(req.sequence, mode=req.mode)
    amps = designer.design()

    # scoring
    for amp in amps:
        score_amplicon(amp, req.blast_db)

    # score 기준 sorting (낮을수록 좋은 점수)
    amps.sort(key=lambda x: x.final_score)

    # 파일 저장
    job_id = uuid.uuid4().hex
    out_json = RESULTS / f"{job_id}.json"

    export = []
    for a in amps:
        export.append({
            "forward": a.forward_primer.sequence,
            "reverse": a.reverse_primer.sequence,
            "probe": a.probe.sequence if a.probe else None,
            "score": a.final_score,
            "f_self_dimer": a.f_self_dimer,
            "r_self_dimer": a.r_self_dimer,
            "hetero_dimer": a.hetero_dimer,
            "f_offtarget": a.f_offtarget,
            "r_offtarget": a.r_offtarget,
            "probe_offtarget": a.probe_offtarget,
        })

    out_json.write_text(json.dumps(export, indent=2))

    return JSONResponse({"job_id": job_id, "result": export})


############################################################
# JSON 다운로드
############################################################

@app.get("/download/json/{job_id}")
def download_json(job_id: str):
    path = RESULTS / f"{job_id}.json"
    if not path.exists():
        raise HTTP