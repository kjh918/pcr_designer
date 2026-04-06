/**
 * qpcr_app.js (Complete Version)
 * - LocalStorage 자동 저장/복구
 * - FastAPI 통신 (최신 Pydantic API Schema 완벽 호환!)
 * - 422 에러 객체 파싱 (object Object 출력 방지)
 * - 26-column Detailed Analytics Table 적용 완료
 */

const state = {
  lastResults: null,
  lastPayload: null,
};

const STORAGE_PREFIX = "qpcr_";

const TRACKED_INPUTS = [
  "input-design-name", "manual-sequence", "input-ref-genome",
  "min_amplicon_length", "max_amplicon_length",
  "primer_min_length", "primer_opt_length", "primer_max_length",
  "primer_min_tm", "primer_opt_tm", "primer_max_tm",
  "primer_min_gc", "primer_opt_gc", "primer_max_gc",
  "probe_min_length", "probe_opt_length", "probe_max_length",
  "probe_min_tm", "probe_opt_tm", "probe_max_tm",
  "probe_min_gc", "probe_opt_gc", "probe_max_gc",
  "qc_hairpin_min_dg", "qc_homodimer_min_dg", "qc_heterodimer_min_dg",
  "qc_min_identity", "qc_min_hit_length", "qc_blast_max_alignments",
  "qc_min_amp_size", "qc_max_amp_size", "qc_use_ispcr_check", 
  "qc_primer_max_diff_tm",
  "qc_probe_poly_g", "qc_probe_avoid_5g"
];

const API_URL = "http://192.168.0.35:9000/api/design/qpcr";

const $ = (sel) => document.querySelector(sel);
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
  if (type === "bool") return raw === "true";

  return raw;
}

function setText(id, text) {
  const el = byId(id);
  if (el) el.innerText = text;
}

function setStatus(text, color) {
  const el = byId("summary-status");
  if (!el) return;
  el.innerText = text;
  if (color) {
      if (text === "SUCCESS" || text === "COMPLETED") {
          el.style.backgroundColor = "#d1e7dd";
          el.style.color = "#0f5132";
      } else if (text.includes("FAIL") || text.includes("ERROR")) {
          el.style.backgroundColor = "#f8d7da";
          el.style.color = "#842029";
      } else {
          el.style.backgroundColor = "#e9ecef";
          el.style.color = "#333";
      }
  }
}

function saveInputsToStorage() {
  TRACKED_INPUTS.forEach((id) => {
    const el = byId(id);
    if (!el) return;
    const val = (el.type === "checkbox") ? el.checked : el.value;
    localStorage.setItem(STORAGE_PREFIX + id, String(val));
  });
}

function restoreInputsFromStorage() {
  TRACKED_INPUTS.forEach((id) => {
    const saved = localStorage.getItem(STORAGE_PREFIX + id);
    const el = byId(id);
    if (!el || saved === null) return;

    if (el.type === "checkbox") el.checked = (saved === "true");
    else el.value = saved;
  });
}

// 🔥 [핵심 보정] QPCRDesignApiInput 모델 구조와 100% 일치하는 계층형 Payload
function buildPayload() {
  return {
    design_name: getVal("input-design-name", "string", "qPCR_Design"),
    sequence: (byId("manual-sequence")?.value || "").replace(/\s/g, "").toUpperCase(),
    reference_genome: getVal("input-ref-genome", "string", "hg38"),
    top_k: 10,
    
    amplicon: {
      min_length: getVal("min_amplicon_length", "int", 60),
      max_length: getVal("max_amplicon_length", "int", 150)
    },
    
    primer: {
      min_length: getVal("primer_min_length", "int", 20),
      opt_length: getVal("primer_opt_length", "int", 25),
      max_length: getVal("primer_max_length", "int", 30),
      min_tm: getVal("primer_min_tm", "float", 55.0),
      opt_tm: getVal("primer_opt_tm", "float", 60.0),
      max_tm: getVal("primer_max_tm", "float", 65.0),
      min_gc: getVal("primer_min_gc", "float", 35.0),
      opt_gc: getVal("primer_opt_gc", "float", 50.0),
      max_gc: getVal("primer_max_gc", "float", 65.0)
    },
    
    probe: {
      min_length: getVal("probe_min_length", "int", 20),
      opt_length: getVal("probe_opt_length", "int", 25),
      max_length: getVal("probe_max_length", "int", 30),
      min_tm: getVal("probe_min_tm", "float", 65.0),
      opt_tm: getVal("probe_opt_tm", "float", 67.0),
      max_tm: getVal("probe_max_tm", "float", 70.0),
      min_gc: getVal("probe_min_gc", "float", 35.0),
      opt_gc: getVal("probe_opt_gc", "float", 50.0),
      max_gc: getVal("probe_max_gc", "float", 65.0),
      max_poly_g: getVal("qc_probe_poly_g", "int", 3),
      max_3_end_gc: getVal("probe_max_3_end_gc", "int", 2),
      avoid_5_prime_g: getVal("qc_probe_avoid_5g", "bool", true)
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
        max_size: getVal("qc_max_amp_size", "int", 300),
        use_ispcr: getVal("qc_use_ispcr_check", "bool", false)
      },
      oligo: {
        max_tm_diff: getVal("qc_primer_max_diff_tm", "float", 3.0)
      }
    }
  };
}

function validatePayload(p) {
  if (!p.sequence || p.sequence.length < 40) {
    return "Please enter a valid DNA sequence (at least 40bp).";
  }
  if (!p.sequence.includes("[") || !p.sequence.includes("]")) {
    return "서열에 대괄호('[', ']')를 사용하여 증폭 타겟을 명시해야 합니다. (e.g. ATGC[A/G]ATGC)";
  }
  return null;
}

function renderSummary(results, payload) {
  const summary = results.summary || {};
  
  setText("summary-project", payload.design_name || "qPCR_Design");
  setText("summary-reference", payload.reference_genome || "hg38");
  
  const now = new Date();
  setText("summary-date", now.toLocaleDateString() + " " + now.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" }));

  const passedCount = summary.passed_count || 0;
  const totalCount = summary.total_count || 0;
  setText("summary-counts", `${passedCount} / ${totalCount}`);

  if (results.status === "success") {
      const allFailed = passedCount === 0 && totalCount > 0;
      setStatus(allFailed ? `FAIL (${passedCount} PASS)` : "SUCCESS", "color");
  } else {
      setStatus("ERROR", "color");
  }

  // 매핑된 계층 구조로 화면 업데이트
  setText("summary-param-amp", `${payload.amplicon.min_length}~${payload.amplicon.max_length} bp`);
  setText("summary-param-primer-tm", `${payload.primer.min_tm}~${payload.primer.max_tm} ℃`);
  setText("summary-param-probe-tm", `${payload.probe.min_tm}~${payload.probe.max_tm} ℃`);
  setText("summary-param-primer-gc", `${payload.primer.min_gc}~${payload.primer.max_gc} %`);

  setText("summary-qc-hairpin", `${payload.qc_criteria.thermodynamics.hairpin_min_dg} kcal`);
  setText("summary-qc-homo", `${payload.qc_criteria.thermodynamics.homodimer_min_dg} kcal`);
  setText("summary-qc-hetero", `${payload.qc_criteria.thermodynamics.heterodimer_min_dg} kcal`);
  setText("summary-qc-tmdiff", `${payload.qc_criteria.oligo.max_tm_diff} ℃`);

  setText("summary-qc-ident", `${payload.qc_criteria.blast.min_identity} %`);
  setText("summary-qc-hitlen", `${payload.qc_criteria.blast.min_hit_length} bp`);
  setText("summary-qc-aligns", `${payload.qc_criteria.blast.max_alignments}`);
  setText("summary-qc-ampsize", `${payload.qc_criteria.amplicon.min_size}~${payload.qc_criteria.amplicon.max_size} bp`);

  setText("summary-qc-probe-5g", payload.probe.avoid_5_prime_g ? "Avoid" : "Allow");
  setText("summary-qc-probe-polyg", `≤ ${payload.probe.max_poly_g}`);
  setText("summary-qc-probe-tmdiff", `Absolute Tm Rules Applied`);
}

function renderResults(results) {
  const tbody = byId("candidate-tbody");
  if (!tbody) return;
  tbody.innerHTML = "";

  const listData = results?.results || [];
  if (!listData.length) {
      tbody.innerHTML = '<tr><td colspan="26" style="text-align:center; padding:60px; color:#888; font-size: 13px;">No candidates found.</td></tr>';
      return;
  }

  listData.forEach((item, idx) => {
      const tr = document.createElement("tr");
      tr.style.cursor = "pointer";
      tr.id = `rank-row-${idx}`;

      const isPass = item.qc_info?.is_pass;
      const qcText = isPass ? "PASS" : (item.qc_info?.fail_reason || "FAIL");
      const qcColor = isPass ? "#198754" : "#d9534f";

      const oligos = item.oligos || {};
      const fwd = oligos.forward || { sequence: "-", tm: 0, gc: 0, cpg_count: 0, hairpin_dg: 0, homodimer_dg: 0 };
      const rev = oligos.reverse || { sequence: "-", tm: 0, gc: 0, cpg_count: 0, hairpin_dg: 0, homodimer_dg: 0 };
      const prb = oligos.probe || { sequence: "-", tm: 0, gc: 0, cpg_count: 0, hairpin_dg: 0, homodimer_dg: 0 };
      const hetero = oligos.heterodimer || { fr_dg: 0, fp_dg: 0, rp_dg: 0 };
      
      const amp = item.amplicon_info || {};
      const size = amp.size ?? "-";
      const ampTm = amp.tm ?? 0.0;
      const pos = amp.genomic_pos ?? "-";
      const alignText = amp.alignment_text_block ?? "No alignment data.";

      tr.innerHTML = `
          <td class="text-center" style="border-right: 1px solid #dee2e6;">${item.rank ?? (idx + 1)}</td>
          <td class="mono-cell">${fwd.sequence}</td>
          <td class="mono-cell">${rev.sequence}</td>
          <td class="mono-cell" style="border-right: 2px solid #adb5bd;">${prb.sequence !== "-" ? prb.sequence : "-"}</td>
          
          <td class="text-center">${formatNum1(fwd.tm)}</td>
          <td class="text-center">${formatNum1(fwd.gc)}%</td>
          <td class="text-center">${fwd.cpg_count ?? 0}</td>
          <td class="text-center" style="color: #d9534f;">${formatNum1(fwd.hairpin_dg)}</td>
          <td class="text-center" style="color: #d9534f; border-right: 2px solid #adb5bd;">${formatNum1(fwd.homodimer_dg)}</td>

          <td class="text-center">${formatNum1(rev.tm)}</td>
          <td class="text-center">${formatNum1(rev.gc)}%</td>
          <td class="text-center">${rev.cpg_count ?? 0}</td>
          <td class="text-center" style="color: #d9534f;">${formatNum1(rev.hairpin_dg)}</td>
          <td class="text-center" style="color: #d9534f; border-right: 2px solid #adb5bd;">${formatNum1(rev.homodimer_dg)}</td>

          <td class="text-center">${formatNum1(prb.tm)}</td>
          <td class="text-center">${formatNum1(prb.gc)}%</td>
          <td class="text-center">${prb.cpg_count ?? 0}</td>
          <td class="text-center" style="color: #d9534f;">${formatNum1(prb.hairpin_dg)}</td>
          <td class="text-center" style="color: #d9534f; border-right: 2px solid #adb5bd;">${formatNum1(prb.homodimer_dg)}</td>

          <td class="text-center" style="color: #d9534f;">${formatNum1(hetero.fr_dg)}</td>
          <td class="text-center" style="color: #d9534f;">${formatNum1(hetero.fp_dg)}</td>
          <td class="text-center" style="color: #d9534f; border-right: 2px solid #adb5bd;">${formatNum1(hetero.rp_dg)}</td>

          <td class="text-center font-weight-bold">${size}</td>
          <td class="text-center font-weight-bold" style="color: #1e707a;">${formatNum1(ampTm)}</td>
          <td class="text-nowrap text-center">${pos}</td>
          <td style="color: ${qcColor}; font-weight: normal; font-size: 11px; white-space: normal; word-break: break-word; line-height: 1.3;">
              ${qcText}
          </td>
      `;

      tr.addEventListener("click", () => {
          document.querySelectorAll("#candidate-tbody tr").forEach((r) => r.classList.remove("selected-row"));
          tr.classList.add("selected-row");
          const alnView = byId("alignment-view");
          if (alnView) alnView.innerText = alignText;
      });

      tbody.appendChild(tr);
  });

  setTimeout(() => {
      const firstRow = document.getElementById("rank-row-0");
      if (firstRow) firstRow.click();
  }, 50);
}

function exportJSON(e) {
  if (e) e.preventDefault();
  if (!state.lastResults) {
    alert("내보낼 데이터가 없습니다. 먼저 분석을 실행하세요.");
    return;
  }

  const projectOutput = state.lastResults;
  const projectName = projectOutput.metadata?.project_name || "qPCR_Project";

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

async function runDesign() {
  const runBtn = byId("btn-run-qpcr") || byId("btn-run-mspcr"); // 범용 실행 버튼
  if (!runBtn) return;

  saveInputsToStorage();
  const payload = buildPayload();
  const err = validatePayload(payload);
  if (err) {
    alert(err);
    return;
  }

  const originalText = runBtn.innerText;
  runBtn.innerText = "⏳ RUNNING...";
  runBtn.disabled = true;
  setStatus("RUNNING", "#2b384c");

  const errPanel = byId("error-log-panel");
  if (errPanel) errPanel.style.display = "none";

  try {
    const res = await fetch(API_URL, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });

    // 🔥 [오브젝트 오브젝트] 출력 방지를 위한 강력한 에러 해석기
    if (!res.ok) {
      const errData = await res.json().catch(() => ({}));
      let errMsg = errData.detail || "Server Error";
      
      if (Array.isArray(errMsg)) {
        errMsg = errMsg.map(e => `[${e.loc?.join('.')}] : ${e.msg}`).join('\n');
      } else if (typeof errMsg === 'object') {
        errMsg = JSON.stringify(errMsg, null, 2);
      }
      throw new Error(errMsg);
    }

    const results = await res.json();
    state.lastPayload = payload;
    state.lastResults = results;

    if (results.status !== "success" || results.error || results.reason) {
        const errContent = byId("error-log-content");
        if (errPanel && errContent) {
            errPanel.style.display = "block";
            let logs = results.log_messages || results.error || results.reason || "Unknown error occurred.";
            if (Array.isArray(logs)) logs = logs.join('\n');
            else if (typeof logs === 'object') logs = JSON.stringify(logs, null, 2);
            errContent.innerText = logs;
        }
    }

    renderResults(results);
    renderSummary(results, payload);

  } catch (e) {
    // 🔥 에러를 깔끔한 텍스트로 알림
    alert("Design failed:\n" + e.message);
    setStatus("ERROR", "red");
    const errContent = byId("error-log-content");
    if (errPanel && errContent) {
        errPanel.style.display = "block";
        errContent.innerText = e.message;
    }
  } finally {
    runBtn.innerText = originalText;
    runBtn.disabled = false;
  }
}

document.addEventListener("DOMContentLoaded", () => {
  restoreInputsFromStorage();

  const runBtn = byId("btn-run-qpcr") || byId("btn-run-mspcr");
  runBtn?.addEventListener("click", (e) => {
    e.preventDefault();
    runDesign();
  });

  byId("btn-export-json")?.addEventListener("click", exportJSON);
  byId("btn-export-html")?.addEventListener("click", exportHTML);

  byId("btn-load-tutorial")?.addEventListener("click", (e) => {
      e.preventDefault();
      if (byId("input-design-name")) byId("input-design-name").value = "EGFR_L858R_Tutorial";
      if (byId("manual-sequence")) byId("manual-sequence").value = "ATGCGTACGTACGTAGCTAGCTAGCATCGATCG[A/G]TACGTAGCTAGCTAGCTAGCATCGATCGA";
      if (byId("input-ref-genome")) byId("input-ref-genome").value = "none";
      
      const btn = e.target;
      btn.innerText = "✅ Loaded!";
      btn.style.backgroundColor = "#28a745";
      setTimeout(() => {
          btn.innerText = "🧪 Tutorial";
          btn.style.backgroundColor = "#17a2b8";
      }, 2000);
  });
});