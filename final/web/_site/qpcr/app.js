/**
 * qpcr_app.js (Complete Version)
 * - LocalStorage 자동 저장/복구
 * - FastAPI 통신 (최적화된 BaseDesignOutput 및 Pydantic 호환 Payload)
 * - 26-column Detailed Analytics Table 적용 완료
 * - Export: JSON / HTML
 */

const state = {
  lastResults: null,   // 서버 응답 전체
  lastPayload: null,   // 마지막 실행 payload
};

const STORAGE_PREFIX = "qpcr_";

const TRACKED_INPUTS = [
  // Basic
  "input-design-name", "manual-sequence", "input-ref-genome",

  // Amplicon
  "min_amplicon_length", "max_amplicon_length",

  // Primer
  "primer_min_length", "primer_opt_length", "primer_max_length",
  "primer_min_tm", "primer_opt_tm", "primer_max_tm",
  "primer_min_gc", "primer_opt_gc", "primer_max_gc",

  // Probe
  "probe_min_length", "probe_opt_length", "probe_max_length",
  "probe_min_tm", "probe_opt_tm", "probe_max_tm",
  "probe_min_gc", "probe_opt_gc", "probe_max_gc",

  // Constraints
  "probe_max_poly_g", "probe_max_3_end_gc", "probe_avoid_5_prime_g",

  // QC (Thermo, Blast, Amp)
  "qc_hairpin_min_dg", "qc_homodimer_min_dg", "qc_heterodimer_min_dg",
  "qc_min_identity", "qc_min_hit_length", "qc_blast_max_alignments",
  "qc_min_amp_size", "qc_max_amp_size", "qc_use_ispcr_check", 
  
  // QC (Primer, Probe diff)
  "qc_primer_min_diff_tm", "qc_primer_max_diff_tm",
  "qc_probe_min_tm_diff", "qc_probe_max_tm_diff", "qc_probe_poly_g", "qc_probe_avoid_5g"
];

// API endpoint (포트와 주소는 환경에 맞게 조정하세요)
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

  if (el.type === "checkbox") {
    const checked = el.checked;
    if (type === "bool") return !!checked;
    return checked;
  }

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

/* -----------------------------
   2) LocalStorage Save/Restore
------------------------------ */

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

/* -----------------------------
   3) Payload Builder + Validation
------------------------------ */

function buildPayload() {
  return {
    design_name: getVal("input-design-name", "string", "qPCR_Design"),
    sequence: (byId("manual-sequence")?.value || "").replace(/\s/g, "").toUpperCase(),
    reference_genome: getVal("input-ref-genome", "string", "hg38"),
    top_k: 10,
    
    // API Schema 요구사항에 맞춘 동적 파라미터 (Dict 주입)
    pcr_params: {
      primer_kwargs: {
        min_amplicon_length: getVal("min_amplicon_length", "int", 60),
        max_amplicon_length: getVal("max_amplicon_length", "int", 150),
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
      probe_kwargs: {
        min_length: getVal("probe_min_length", "int", 20),
        opt_length: getVal("probe_opt_length", "int", 25),
        max_length: getVal("probe_max_length", "int", 30),
        min_tm: getVal("probe_min_tm", "float", 65.0),
        opt_tm: getVal("probe_opt_tm", "float", 67.0),
        max_tm: getVal("probe_max_tm", "float", 70.0),
        min_gc: getVal("probe_min_gc", "float", 35.0),
        opt_gc: getVal("probe_opt_gc", "float", 50.0),
        max_gc: getVal("probe_max_gc", "float", 65.0),
        max_probe_poly_g: getVal("probe_max_poly_g", "int", 3),
        max_probe_3_end_gc: getVal("probe_max_3_end_gc", "int", 2),
        avoid_5_prime_g: getVal("probe_avoid_5_prime_g", "bool", true)
      }
    },
    
    qc_criteria: {
      hairpin_min_dg: getVal("qc_hairpin_min_dg", "float", -5.0),
      homodimer_min_dg: getVal("qc_homodimer_min_dg", "float", -6.0),
      heterodimer_min_dg: getVal("qc_heterodimer_min_dg", "float", -6.0),
      min_identity: getVal("qc_min_identity", "float", 90.0),
      min_hit_length: getVal("qc_min_hit_length", "int", 13),
      blast_max_alignments: getVal("qc_blast_max_alignments", "int", 50),
      min_amp_size: getVal("qc_min_amp_size", "int", 50),
      max_amp_size: getVal("qc_max_amp_size", "int", 300),
      use_ispcr_check: getVal("qc_use_ispcr_check", "bool", false),
      primer: {
          min_diff_tm: getVal("qc_primer_min_diff_tm", "float", 0.0),
          max_diff_tm: getVal("qc_primer_max_diff_tm", "float", 3.0)
      },
      probe: {
          min_primer_probe_tm_diff: getVal("qc_probe_min_tm_diff", "float", 5.0),
          max_primer_probe_tm_diff: getVal("qc_probe_max_tm_diff", "float", 10.0),
          max_probe_poly_g: getVal("qc_probe_poly_g", "int", 3),
          avoid_5_prime_g: getVal("qc_probe_avoid_5g", "bool", true)
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

/* -----------------------------
   4) Render: Results / Summary
------------------------------ */

function renderSummary(results, payload) {
  const meta = results.metadata || results.export_meta || {};
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

  // 1. Design Target & Params
  setText("summary-param-amp", `${payload.pcr_params.primer_kwargs.min_amplicon_length}~${payload.pcr_params.primer_kwargs.max_amplicon_length} bp`);
  setText("summary-param-primer-tm", `${payload.pcr_params.primer_kwargs.min_tm}~${payload.pcr_params.primer_kwargs.max_tm} ℃`);
  setText("summary-param-probe-tm", `${payload.pcr_params.probe_kwargs.min_tm}~${payload.pcr_params.probe_kwargs.max_tm} ℃`);
  setText("summary-param-primer-gc", `${payload.pcr_params.primer_kwargs.min_gc}~${payload.pcr_params.primer_kwargs.max_gc} %`);

  // 2. Applied QC Parameters
  setText("summary-qc-hairpin", `${payload.qc_criteria.hairpin_min_dg} kcal`);
  setText("summary-qc-homo", `${payload.qc_criteria.homodimer_min_dg} kcal`);
  setText("summary-qc-hetero", `${payload.qc_criteria.heterodimer_min_dg} kcal`);
  setText("summary-qc-tmdiff", `${payload.qc_criteria.primer.max_diff_tm} ℃`);

  setText("summary-qc-ident", `${payload.qc_criteria.min_identity} %`);
  setText("summary-qc-hitlen", `${payload.qc_criteria.min_hit_length} bp`);
  setText("summary-qc-aligns", `${payload.qc_criteria.blast_max_alignments}`);
  setText("summary-qc-ampsize", `${payload.qc_criteria.min_amp_size}~${payload.qc_criteria.max_amp_size} bp`);

  setText("summary-qc-probe-5g", payload.qc_criteria.probe.avoid_5_prime_g ? "Avoid" : "Allow");
  setText("summary-qc-probe-polyg", `≤ ${payload.qc_criteria.probe.max_probe_poly_g}`);
  setText("summary-qc-probe-tmdiff", `${payload.qc_criteria.probe.min_primer_probe_tm_diff}~${payload.qc_criteria.probe.max_primer_probe_tm_diff} ℃`);
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

      // 26-column HTML 렌더링
      tr.innerHTML = `
          <!-- Basic Information -->
          <td class="text-center" style="border-right: 1px solid #dee2e6;">${item.rank ?? (idx + 1)}</td>
          <td class="mono-cell">${fwd.sequence}</td>
          <td class="mono-cell">${rev.sequence}</td>
          <td class="mono-cell" style="border-right: 2px solid #adb5bd;">${prb.sequence !== "-" ? prb.sequence : "-"}</td>
          
          <!-- Forward Primer Details -->
          <td class="text-center">${formatNum1(fwd.tm)}</td>
          <td class="text-center">${formatNum1(fwd.gc)}%</td>
          <td class="text-center">${fwd.cpg_count ?? 0}</td>
          <td class="text-center" style="color: #d9534f;">${formatNum1(fwd.hairpin_dg)}</td>
          <td class="text-center" style="color: #d9534f; border-right: 2px solid #adb5bd;">${formatNum1(fwd.homodimer_dg)}</td>

          <!-- Reverse Primer Details -->
          <td class="text-center">${formatNum1(rev.tm)}</td>
          <td class="text-center">${formatNum1(rev.gc)}%</td>
          <td class="text-center">${rev.cpg_count ?? 0}</td>
          <td class="text-center" style="color: #d9534f;">${formatNum1(rev.hairpin_dg)}</td>
          <td class="text-center" style="color: #d9534f; border-right: 2px solid #adb5bd;">${formatNum1(rev.homodimer_dg)}</td>

          <!-- Probe Details -->
          <td class="text-center">${formatNum1(prb.tm)}</td>
          <td class="text-center">${formatNum1(prb.gc)}%</td>
          <td class="text-center">${prb.cpg_count ?? 0}</td>
          <td class="text-center" style="color: #d9534f;">${formatNum1(prb.hairpin_dg)}</td>
          <td class="text-center" style="color: #d9534f; border-right: 2px solid #adb5bd;">${formatNum1(prb.homodimer_dg)}</td>

          <!-- Heterodimer (dG) -->
          <td class="text-center" style="color: #d9534f;">${formatNum1(hetero.fr_dg)}</td>
          <td class="text-center" style="color: #d9534f;">${formatNum1(hetero.fp_dg)}</td>
          <td class="text-center" style="color: #d9534f; border-right: 2px solid #adb5bd;">${formatNum1(hetero.rp_dg)}</td>

          <!-- Amplicon & QC -->
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

  // 첫 행 자동 선택
  setTimeout(() => {
      const firstRow = document.getElementById("rank-row-0");
      if (firstRow) firstRow.click();
  }, 50);
}

/* -----------------------------
   5) Export: JSON / HTML
------------------------------ */

function exportJSON(e) {
  if (e) e.preventDefault();
  if (!state.lastResults) {
    alert("내보낼 데이터가 없습니다. 먼저 분석을 실행하세요.");
    return;
  }

  const projectOutput = state.lastResults;
  const projectName = projectOutput.export_meta?.project_name || "qPCR_Project";

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

/* -----------------------------
   6) Main Run (Fetch)
------------------------------ */

async function runDesign() {
  const runBtn = byId("btn-run-qpcr"); // 사용하시는 HTML의 실행버튼 ID에 맞추세요
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

    if (!res.ok) {
      const errData = await res.json().catch(() => ({}));
      throw new Error(errData.detail || "Server Error");
    }

    const results = await res.json();
    state.lastPayload = payload;
    state.lastResults = results;

    if (results.status !== "success" || results.error || results.reason) {
        const errContent = byId("error-log-content");
        if (errPanel && errContent) {
            errPanel.style.display = "block";
            errContent.innerText = results.log_messages || results.error || results.reason || "Unknown error occurred.";
        }
    }

    renderResults(results);
    renderSummary(results, payload);

  } catch (e) {
    alert("Design failed: " + e.message);
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

/* -----------------------------
   7) Wire up Events
------------------------------ */

document.addEventListener("DOMContentLoaded", () => {
  restoreInputsFromStorage();

  byId("btn-run-qpcr")?.addEventListener("click", (e) => {
    e.preventDefault();
    runDesign();
  });

  byId("btn-export-json")?.addEventListener("click", exportJSON);
  byId("btn-export-html")?.addEventListener("click", exportHTML);

  // Tutorial 셋업 (예: EGFR)
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