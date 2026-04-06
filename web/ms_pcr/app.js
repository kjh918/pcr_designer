/**
 * [v5.0] MS-PCR app.js (Integrated Architecture)
 * - AS-PCR 스타일의 Sets > Alleles (M/U) 구조 완벽 지원
 * - 행 클릭 시 alignment_text_block 실시간 바인딩
 * - 요약 패널(Summary) 및 422 에러 상세 파싱 포함
 */

const state = {
    lastResults: null,
    lastPayload: null,
};

const STORAGE_PREFIX = "mspcr_";
const API_URL = `http://${window.location.hostname}:9000/api/design/mspcr`;

const byId = (id) => document.getElementById(id);

// --- 헬퍼 함수 ---
function setText(id, text) {
    const el = byId(id);
    if (el) el.innerText = text;
}

function formatNum1(val) {
    if (val === "-" || val === undefined || val === null || Number.isNaN(val)) return "-";
    const n = typeof val === "number" ? val : parseFloat(val);
    return Number.isFinite(n) ? n.toFixed(1) : "-";
}

function getVal(id, type = "string", def = null) {
    const el = byId(id);
    if (!el) return def;
    if (el.type === "checkbox") return type === "bool" ? !!el.checked : el.checked;
    const raw = (el.value ?? "").toString();
    if (type === "int") return parseInt(raw, 10) || 0;
    if (type === "float") return parseFloat(raw) || 0.0;
    return raw;
}

function setStatus(text) {
    const el = byId("summary-status");
    if (!el) return;
    el.innerText = text;
    if (text === "SUCCESS") { el.style.backgroundColor = "#d1e7dd"; el.style.color = "#0f5132"; }
    else if (text === "RUNNING") { el.style.backgroundColor = "#fff3cd"; el.style.color = "#856404"; }
    else if (text === "ERROR") { el.style.backgroundColor = "#f8d7da"; el.style.color = "#842029"; }
}

// --- 1. 디자인 요약 패널 업데이트 (AS-PCR 스타일) ---
function renderSummary(results, payload) {
    const metadata = results.metadata || {};
    const sets = results.results || [];
    
    setText("summary-project", payload.design_name || "MSPCR_Project");
    setText("summary-reference", metadata.reference_genome || payload.reference_genome || "hg38");
    setText("summary-date", new Date().toLocaleDateString() + " " + new Date().toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" }));

    if (sets.length > 0) {
        const firstSet = sets[0];
        const anyAllele = firstSet.alleles?.M || firstSet.alleles?.U;
        if (anyAllele) {
            setText("summary-pos", anyAllele.amplicon_info?.genomic_pos || "-");
            const pos = anyAllele.amplicon_info?.genomic_pos || "";
            if (pos.includes(":")) setText("summary-chrom", pos.split(":")[0]);
        }
        setText("summary-mutation", "[CG] Methylation Target");
    }
}

// --- 2. MS-PCR 전용 Payload 빌더 (계층형 구조) ---
function buildPayload() {
    const rawSeq = (byId("mspcr-raw-sequence")?.value || "").replace(/\s/g, "").toUpperCase();
    return {
        design_name: getVal("mspcr-design-name") || "MSPCR_Project",
        sequence: rawSeq,
        reference_genome: getVal("input-ref-genome", "string", "hg38"),
        top_k: 5,
        window_size_3prime: 3, // MS-PCR 특화
        min_cpg_count: 1,      // MS-PCR 특화
        amplicon: {
            min_length: getVal("min_amplicon_length", "int", 80),
            max_length: getVal("max_amplicon_length", "int", 200)
        },
        primer: {
            min_length: getVal("primer_min_length", "int", 15),
            opt_length: getVal("primer_opt_length", "int", 18),
            max_length: getVal("primer_max_length", "int", 25),
            min_tm: getVal("primer_min_tm", "float", 50.0),
            opt_tm: getVal("primer_opt_tm", "float", 55.0),
            max_tm: getVal("primer_max_tm", "float", 60.0),
            min_gc: getVal("primer_min_gc", "float", 10.0),
            opt_gc: getVal("primer_opt_gc", "float", 25.0),
            max_gc: getVal("primer_max_gc", "float", 60.0)
        },
        qc_criteria: {
            thermodynamics: {
                hairpin_min_dg: getVal("qc_hairpin_min_dg", "float", -5.0),
                homodimer_min_dg: getVal("qc_homodimer_min_dg", "float", -6.0),
                heterodimer_min_dg: getVal("qc_heterodimer_min_dg", "float", -6.0)
            },
            blast: {
                min_identity: getVal("qc_min_identity", "float", 90.0),
                min_hit_length: getVal("qc_min_hit_length", "int", 13),
                max_alignments: getVal("qc_blast_max_alignments", "int", 50)
            },
            amplicon: {
                min_size: getVal("qc_min_amp_size", "int", 50),
                max_size: getVal("qc_max_amp_size", "int", 300)
            },
            oligo: {
                max_tm_diff: getVal("qc_primer_max_diff_tm", "float", 3.0)
            }
        }
    };
}

// --- 3. 정렬 정보 실시간 바인딩 (AS-PCR 스타일) ---
window.showMspcrAlignment = function(alleleId) {
    if (!state.lastResults || !state.lastResults.results) return;
    let foundAllele = null;
    state.lastResults.results.forEach(set => {
        if (set.alleles) {
            Object.values(set.alleles).forEach(allele => {
                if (allele.id === alleleId) foundAllele = allele;
            });
        }
    });
    const view = byId("alignment-view");
    if (view && foundAllele?.amplicon_info?.alignment_text_block) {
        view.textContent = foundAllele.amplicon_info.alignment_text_block;
        document.querySelectorAll(".allele-detail-row").forEach(r => r.classList.remove("selected-row"));
        const targetRow = document.querySelector(`tr[data-allele-id="${alleleId}"]`);
        if (targetRow) targetRow.classList.add("selected-row");
    }
};

function exportJSON(e) {
  if (e) e.preventDefault();
  if (!state.lastResults) {
    alert("내보낼 데이터가 없습니다. 먼저 분석을 실행하세요.");
    return;
  }

  const projectOutput = state.lastResults;
  const projectName = projectOutput.metadata?.project_name || "MSPCR_Project";

  const jsonString = JSON.stringify(projectOutput, null, 4);
  const blob = new Blob([jsonString], { type: "application/json" });
  const url = URL.createObjectURL(blob);

  const link = document.createElement("a");
  link.href = url;
  link.download = `${projectName}_QC_Report.json`;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}

function exportHTML(e) {
  if (e) e.preventDefault();
  if (!state.lastResults) {
      alert("내보낼 데이터가 없습니다. 먼저 분석을 실행하세요.");
      return;
  }

  const target = document.querySelector(".main-content");
  if (!target) {
      alert("리포트 영역을 찾을 수 없습니다.");
      return;
  }

  let styles = "";
  try {
      styles = Array.from(document.styleSheets)
          .map((ss) => {
              try { return Array.from(ss.cssRules).map((r) => r.cssText).join("\n"); } 
              catch { return ""; }
          }).join("\n");
  } catch(err) {
      console.warn("CSS styles could not be fully exported.", err);
  }

  const projectName = state.lastPayload?.design_name || "qPCR_Project";
  const htmlContent = `<!DOCTYPE html>
<html lang="ko">
<head>
  <meta charset="UTF-8" />
  <title>qPCR Design Report - ${projectName}</title>
  <style>
    body { font-family: 'Nunito', sans-serif; padding: 20px; background: #f5f5f5; }
    ${styles}
    .no-print, details { display: none !important; }
    .panel { background: white; padding: 20px; border-radius: 8px; border: 1px solid #ddd; margin-bottom: 20px; }
    .table-container { max-height: none !important; overflow: visible !important; }
    pre { background: #2b2b2b; color: #a9b7c6; padding: 15px; border-radius: 4px; }
  </style>
</head>
<body>
  <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1);">
    <h2 style="color: #1e707a; margin-top: 0;">🔍 qPCR Design & Analysis Report</h2>
    <p style="color: #666; font-size: 13px;">Export Date: ${new Date().toLocaleString()}</p>
    <hr style="border: none; border-bottom: 2px solid #eee; margin-bottom: 20px;">
    ${target.innerHTML}
  </div>
</body>
</html>`;

  const blob = new Blob([htmlContent], { type: "text/html;charset=utf-8;" });
  const url = URL.createObjectURL(blob);

  const link = document.createElement("a");
  link.href = url;
  link.download = `${projectName}_Report.html`;
  link.click();
  URL.revokeObjectURL(url);
}

// --- 4. 결과 테이블 렌더링 (Sets > M/U Alleles 구조) ---
function renderResults(results) {
    const tbody = byId("candidate-tbody");
    if (!tbody) return;
    tbody.innerHTML = "";

    const sets = results?.results || [];
    if (!sets.length) {
        tbody.innerHTML = '<tr><td colspan="17" style="text-align:center; padding:60px; color:#999;">No candidates found.</td></tr>';
        return;
    }

    sets.forEach((set, setIdx) => {
        ["M", "U"].forEach((key) => {
            const data = set.alleles[key];
            if (!data) return;

            const tr = document.createElement("tr");
            tr.className = "allele-detail-row";
            tr.style.cursor = "pointer";
            tr.setAttribute("data-allele-id", data.id);

            const isPass = data.qc_info?.is_pass;
            const qcColor = isPass ? "#198754" : "#d9534f";
            const labelBg = key === "M" ? "#d9534f" : "#0275d8";
            
            tr.innerHTML = `
                <td class="text-center">${key === "M" ? (set.rank || setIdx + 1) : ""}</td>
                <td class="text-center"><span style="background:${labelBg}; color:white; padding:2px 6px; border-radius:3px; font-size:10px; font-weight:bold;">${key}-Allele</span></td>
                <td class="text-center" style="border-right: 2px solid #adb5bd; font-weight:bold;">${key === "M" ? set.set_id : ""}</td>
                
                <td class="mono-cell">${data.oligos.forward.sequence}</td>
                <td class="text-center">${formatNum1(data.oligos.forward.tm)}</td>
                <td class="text-center">${formatNum1(data.oligos.forward.gc)}</td>
                <td class="text-center" style="color:#d9534f;">${formatNum1(data.oligos.forward.hairpin_dg)}</td>
                <td class="text-center" style="color:#d9534f; border-right: 2px solid #adb5bd;">${formatNum1(data.oligos.forward.homodimer_dg)}</td>

                <td class="mono-cell">${data.oligos.reverse.sequence}</td>
                <td class="text-center">${formatNum1(data.oligos.reverse.tm)}</td>
                <td class="text-center">${formatNum1(data.oligos.reverse.gc)}</td>
                <td class="text-center" style="color:#d9534f;">${formatNum1(data.oligos.reverse.hairpin_dg)}</td>
                <td class="text-center" style="color:#d9534f; border-right: 2px solid #adb5bd;">${formatNum1(data.oligos.reverse.homodimer_dg)}</td>

                <td class="text-center" style="color:#d9534f; border-right: 2px solid #adb5bd;">${formatNum1(data.oligos.heterodimer?.fr_dg || 0)}</td>
                <td class="text-center">${formatNum1(data.metrics.pair_penalty)}</td>
                <td class="text-center" style="font-size: 10px;">${data.amplicon_info.genomic_pos}</td>
                <td style="color: ${qcColor}; font-size: 10px; line-height: 1.2;">${data.qc_info.fail_reason}</td>
            `;

            tr.addEventListener("click", () => window.showMspcrAlignment(data.id));
            tbody.appendChild(tr);
        });
    });

    // 자동 첫 번째 행 선택
    setTimeout(() => {
        const firstRow = tbody.querySelector(".allele-detail-row");
        if (firstRow) firstRow.click();
    }, 150);
}

// --- 5. 분석 실행 메인 ---
async function runDesign() {
    const runBtn = byId("btn-run-mspcr");
    if (!runBtn) return;

    const payload = buildPayload();
    if (!payload.sequence.includes("[") || !payload.sequence.includes("]")) {
        return alert("타겟 CpG를 [CG]로 표시해주세요.");
    }

    runBtn.innerText = "⏳ RUNNING...";
    runBtn.disabled = true;
    setStatus("RUNNING");

    try {
        const res = await fetch(API_URL, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(payload),
        });

        if (!res.ok) {
            const errData = await res.json().catch(() => ({}));
            let msg = errData.detail || "Server Error";
            if (Array.isArray(msg)) msg = msg.map(e => `[${(e.loc || []).join('.')}] ${e.msg}`).join('\n');
            throw new Error(msg);
        }

        const results = await res.json();
        state.lastResults = results;
        state.lastPayload = payload;

        renderResults(results);
        renderSummary(results, payload);
        setStatus("SUCCESS");

    } catch (e) {
        alert("Design failed: " + e.message);
        setStatus("ERROR");
    } finally {
        runBtn.innerText = "DESIGN START";
        runBtn.disabled = false;
    }
}

// --- 6. 초기화 및 이벤트 등록 ---
document.addEventListener("DOMContentLoaded", () => {
    // 저장된 입력 복구
    const inputs = ["mspcr-design-name", "mspcr-raw-sequence", "input-ref-genome"];
    inputs.forEach(id => {
        const saved = localStorage.getItem(STORAGE_PREFIX + id);
        if (saved && byId(id)) byId(id).value = saved;
    });

    byId("btn-export-json")?.addEventListener("click", exportJSON);
    byId("btn-export-html")?.addEventListener("click", exportHTML);
    byId("btn-run-mspcr")?.addEventListener("click", (e) => { e.preventDefault(); runDesign(); });
    
    window.addEventListener("beforeunload", () => {
        inputs.forEach(id => {
            if (byId(id)) localStorage.setItem(STORAGE_PREFIX + id, byId(id).value);
        });
    });
});