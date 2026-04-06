/**
 * app.js (AS-PCR Accordion UI & Unified Params Final Version)
 * - TypeError 방지 로직 및 ID 정합성 강화
 * - API URL 동적 할당 및 422 에러 상세 파싱
 */

const state = {
    lastResults: null,
    lastPayload: null,
};

const STORAGE_PREFIX = "aspcr_";

// 🔥 실제 UI에 존재하는 모든 ID를 빠짐없이 추적합니다.
const TRACKED_INPUTS = [
    "aspcr-design-name", "mspcr-design-name", "input-design-name", 
    "aspcr-raw-sequence", "mspcr-raw-sequence", "manual-sequence", 
    "input-ref-genome",
    "aspcr_fixed_prime", "aspcr_mismatch_pos", "aspcr_mismatch_intensity",
    "min_amplicon_length", "max_amplicon_length",
    "primer_min_length", "primer_opt_length", "primer_max_length",
    "primer_min_tm", "primer_opt_tm", "primer_max_tm",
    "primer_min_gc", "primer_opt_gc", "primer_max_gc",
    "qc_hairpin_min_dg", "qc_homodimer_min_dg", "qc_heterodimer_min_dg",
    "qc_min_identity", "qc_min_hit_length", "qc_blast_max_alignments",
    "qc_min_amp_size", "qc_max_amp_size", "qc_use_ispcr_check", 
    "qc_primer_max_diff_tm"
];

// 🔥 접속한 서버 주소를 자동으로 찾아가도록 유동적 API 설정
const API_URL = `${window.location.protocol}//${window.location.hostname}:9000/api/design/aspcr`;

const byId = (id) => document.getElementById(id);

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
    if (type === "int") {
        const v = parseInt(raw, 10);
        return Number.isFinite(v) ? v : def;
    }
    if (type === "float") {
        const v = parseFloat(raw);
        return Number.isFinite(v) ? v : def;
    }
    return raw;
}
function setText(id, text) {
    const el = byId(id);
    if (el) el.innerText = text;
}
function setStatus(text) {
    const el = byId("summary-status");
    if (!el) return;
    el.innerText = text;
    // 상태에 따른 색상 변경
    if (text === "SUCCESS") { el.style.backgroundColor = "#d1e7dd"; el.style.color = "#0f5132"; }
    else if (text === "RUNNING") { el.style.backgroundColor = "#fff3cd"; el.style.color = "#856404"; }
    else if (text === "ERROR") { el.style.backgroundColor = "#f8d7da"; el.style.color = "#842029"; }
}

function buildPayload() {
    // 🔥 ID 우선순위: aspcr- -> mspcr- -> manual-
    const rawSeq = (byId("aspcr-raw-sequence")?.value || byId("mspcr-raw-sequence")?.value || byId("manual-sequence")?.value || "").replace(/\s/g, "").toUpperCase();
    const designName = getVal("aspcr-design-name") || getVal("mspcr-design-name") || getVal("input-design-name") || "ASPCR_Project";

    return {
        design_name: designName,
        sequence: rawSeq,
        reference_genome: getVal("input-ref-genome", "string", "hg38"),
        top_k: 5,
        fixed_prime: getVal("input-strand"), // 이제 직접 "forward" 또는 "reverse"를 가져옴
        mismatch_pos: parseInt(getVal("input-mismatch")), // "0", "3", "2" 등을 숫자로 변환
        mismatch_intensity: getVal("input-mismatch-intensity"), // 직접 "strong", "weak"를 가져옴
        amplicon: {
            min_length: getVal("min_amplicon_length", "int", 60),
            max_length: getVal("max_amplicon_length", "int", 150)
        },
        primer: {
            min_length: getVal("primer_min_length", "int", 15),
            opt_length: getVal("primer_opt_length", "int", 22),
            max_length: getVal("primer_max_length", "int", 30),
            min_tm: getVal("primer_min_tm", "float", 52.0),
            opt_tm: getVal("primer_opt_tm", "float", 58.0),
            max_tm: getVal("primer_max_tm", "float", 65.0),
            min_gc: getVal("primer_min_gc", "float", 35.0),
            opt_gc: getVal("primer_opt_gc", "float", 50.0),
            max_gc: getVal("primer_max_gc", "float", 65.0)
        },
        qc_criteria: {
            thermodynamics: {
                hairpin_min_dg: getVal("qc_hairpin_min_dg", "float", -6.0),
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
                max_size: getVal("qc_max_amp_size", "int", 300),
                use_ispcr: getVal("qc_use_ispcr_check", "bool", false)
            },
            oligo: {
                max_tm_diff: getVal("qc_primer_max_diff_tm", "float", 5.0)
            }
        }
    };
}

function exportJSON(e) {
  if (e) e.preventDefault();
  if (!state.lastResults) {
    alert("내보낼 데이터가 없습니다. 먼저 분석을 실행하세요.");
    return;
  }

  const projectOutput = state.lastResults;
  const projectName = projectOutput.metadata?.project_name || "ASPCR_Project";

  const jsonString = JSON.stringify(projectOutput, null, 4);
  const blob = new Blob([jsonString], { type: "application/json" });
  const url = URL.createObjectURL(blob);

  const link = document.createElement("a");
  link.href = url;
  link.download = `${projectName}_ASPCR_Report.json`;
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

function validatePayload(p) {
    if (!p.sequence || p.sequence.length < 30) return "Please enter a valid DNA sequence (at least 30bp).";
    if (!p.sequence.includes("[") || !p.sequence.includes("]")) return "서열에 대괄호('[', ']')를 사용하여 변이 위치를 표시하세요. (예: ATGC[G,A]ATGC)";
    return null;
}

window.toggleSet = function(setId) {
    const rows = document.querySelectorAll('.' + setId);
    const icon = document.getElementById('icon-' + setId);
    if (!rows.length) return;
    
    let isNowVisible = false;
    rows.forEach(row => {
        if (row.style.display === 'none') {
            row.style.display = 'table-row';
            isNowVisible = true;
        } else {
            row.style.display = 'none';
        }
    });
    if (icon) icon.textContent = isNowVisible ? '-' : '+';
};

function renderSummary(results, payload) {
    const metadata = results.metadata || {};
    const sets = results.results || [];
    
    setText("summary-project", payload.design_name || "qPCR_Design");
    setText("summary-reference", metadata.reference_genome || payload.reference_genome || "hg38");
    
    const now = new Date();
    setText("summary-date", now.toLocaleDateString() + " " + now.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" }));

    // 3. Target Position & Mutation (첫 번째 결과 기준)
    if (sets.length > 0) {
        const firstSet = sets[0];
        const anyAllele = firstSet.alleles?.wt || firstSet.alleles?.alt;
        
        if (anyAllele) {
            // 좌표 정보 채우기
            setText("summary-pos", anyAllele.amplicon_info?.genomic_pos || "-");
        }
        
        // 4. Mutation 정보 추출 (예: [A,G] -> A > G)
        const rawSeq = payload.sequence || "";
        const match = rawSeq.match(/\[(.*?)\]/);
        if (match) {
            const parts = match[1].split(/[,\/]/);
            const ref = parts[0]?.trim() || "?";
            const alt = parts[1]?.trim() || "?";
            setText("summary-mutation", `${ref} > ${alt}`);
        } else {
            setText("summary-mutation", metadata.mutation || "-");
        }
    }
}
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
        const alleleKeys = ["wt", "alt", "wt_mm", "alt_mm"];
        
        alleleKeys.forEach((key) => {
            const data = set.alleles[key];
            if (!data) return;

            const tr = document.createElement("tr");
            tr.className = "allele-detail-row";
            tr.style.cursor = "pointer";
            tr.setAttribute("data-allele-id", data.id);

            const isPass = data.qc_info?.is_pass;
            const qcText = isPass ? "PASS" : (data.qc_info?.fail_reason || "FAIL");
            const qcColor = isPass ? "#198754" : "#d9534f";
            const labelColor = key.includes('mm') ? '#d9534f' : '#0275d8';
            
            const alignText = data.amplicon_info?.alignment_text_block || "No alignment data available.";

            tr.innerHTML = `
                <td class="text-center">${key === "wt" ? (set.rank || setIdx + 1) : ""}</td>
                <td style="font-weight: bold; color: ${labelColor}">${key.toUpperCase()}</td>
                <td class="text-center" style="border-right: 2px solid #adb5bd; font-weight:bold;">${key === "wt" ? set.set_id : ""}</td>
                
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
                <td style="color: ${qcColor}; font-size: 10px; line-height: 1.2;">${qcText}</td>
            `;

            // 🔥 [클릭 이벤트] 정렬 뷰어 업데이트
            tr.addEventListener("click", () => {
                document.querySelectorAll(".allele-detail-row").forEach(r => r.classList.remove("selected-row"));
                tr.classList.add("selected-row");
                const view = byId("alignment-view");
                if (view) view.textContent = alignText;
            });

            tbody.appendChild(tr);
        });
    });

    // 첫 번째 결과 자동 선택
    setTimeout(() => {
        const firstRow = tbody.querySelector(".allele-detail-row");
        if (firstRow) firstRow.click();
    }, 100);
}
async function runDesign() {
    const runBtn = byId("btn-run-aspcr");
    if (!runBtn) return;

    const payload = buildPayload();
    const err = validatePayload(payload);
    if (err) return alert(err);

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
            let errorMsg = "Server Error";
            if (errData.detail) {
                // 🔥 [TypeError 방지] loc가 존재할 때만 join 수행하는 안전한 파싱
                errorMsg = Array.isArray(errData.detail) 
                    ? errData.detail.map(e => `[${(e.loc || []).join('.')}] ${e.msg}`).join('\n') 
                    : errData.detail;
            }
            throw new Error(errorMsg);
        }

        const results = await res.json();
        state.lastResults = results;
        renderResults(results);
        renderSummary(results,payload)
        setStatus("SUCCESS");
        setTimeout(() => window.toggleSet("set-row-1"), 100);

    } catch (e) {
        alert("Design failed: " + e.message);
        setStatus("ERROR");
    } finally {
        runBtn.innerText = "DESIGN START";
        runBtn.disabled = false;
    }
}

document.addEventListener("DOMContentLoaded", () => {
    // LocalStorage 복구
    TRACKED_INPUTS.forEach(id => {
        const saved = localStorage.getItem(STORAGE_PREFIX + id);
        const el = byId(id);
        if (el && saved !== null) el.type === "checkbox" ? el.checked = (saved === "true") : el.value = saved;
    });

    // 이벤트 리스너 등록 (null 체크 포함)
    byId("btn-export-json")?.addEventListener("click", exportJSON);
    byId("btn-export-html")?.addEventListener("click", exportHTML);
    byId("btn-run-aspcr")?.addEventListener("click", (e) => { e.preventDefault(); runDesign(); });
    
    // 자동 저장
    window.addEventListener("beforeunload", () => {
        TRACKED_INPUTS.forEach(id => {
            const el = byId(id);
            if (el) localStorage.setItem(STORAGE_PREFIX + id, el.type === "checkbox" ? String(el.checked) : el.value);
        });
    });
});