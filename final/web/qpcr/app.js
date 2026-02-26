/**
 * qpcr_app.js (refactor)
 * - LocalStorage 자동 저장/복구
 * - FastAPI 통신
 * - 결과/메타 렌더링
 * - Export: JSON / HTML (PDF는 별도 구현 가능)
 */

/* -----------------------------
   0) App State / Constants
------------------------------ */

const state = {
  lastResults: null,   // 서버 응답 전체
  lastPayload: null,   // 마지막 실행 payload
};

const STORAGE_PREFIX = "qpcr_";

const TRACKED_INPUTS = [
  // Basic
  "input-ref-genome", "input-chrom", "input-start", "input-end",
  "input-ref-base", "input-alt-base", "input-strand",

  // Amplicon
  "min_amplicon_length", "max_amplicon_length",

  // Primer
  "primer_min_length", "primer_opt_length", "primer_max_length",
  "primer_min_tm", "primer_opt_tm", "primer_max_tm",
  "primer_min_gc", "primer_opt_gc", "primer_max_gc",

  // Probe
  "probe_min_length", "probe_opt_length", "probe_max_length",
  "min_primer_probe_tm_diff", "max_primer_probe_tm_diff",
  "probe_min_tm", "probe_opt_tm", "probe_max_tm",
  "probe_min_gc", "probe_opt_gc", "probe_max_gc",

  // Constraints
  "probe_max_poly_g", "probe_max_3_end_gc", "probe_avoid_5_prime_g",

  // QC
  "qc_hairpin_min_dg", "qc_homodimer_min_dg", "qc_heterodimer_min_dg",
  "qc_min_identity", "qc_min_hit_length", "qc_blast_max_alignments",
  "qc_use_ispcr_check", "qc_primer_max_diff_tm",
];

// API endpoint (필요하면 환경에 맞게 바꾸기)
const API_URL = "http://192.168.0.35:9000/api/design/qpcr";

/* -----------------------------
   1) Utils (DOM / format / read)
------------------------------ */

const $ = (sel) => document.querySelector(sel);

function byId(id) {
  return document.getElementById(id);
}

function formatNum1(val) {
  if (val === "-" || val === undefined || val === null || Number.isNaN(val)) return "-";
  const n = typeof val === "number" ? val : parseFloat(val);
  return Number.isFinite(n) ? n.toFixed(1) : "-";
}

function getVal(id, type = "string", def = null) {
  const el = byId(id);
  if (!el) return def;

  // checkbox는 checked 기준
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
  if (type === "bool") {
    // select/hidden 등에서 "true"/"false"로 저장되는 경우
    return raw === "true";
  }

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
  if (color) el.style.color = color;
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
  console.log("💾 Input values saved.");
}

function restoreInputsFromStorage() {
  TRACKED_INPUTS.forEach((id) => {
    const saved = localStorage.getItem(STORAGE_PREFIX + id);
    const el = byId(id);
    if (!el || saved === null) return;

    if (el.type === "checkbox") {
      el.checked = (saved === "true");
    } else {
      el.value = saved;
    }
  });
  console.log("📂 Input values restored.");
}

/* -----------------------------
   3) Payload Builder + Validation
------------------------------ */

function buildPayload() {
  return {
    // Basic
    reference: getVal("input-ref-genome", "string", "hg38"),
    chrom: getVal("input-chrom", "string", "").trim(),
    start: getVal("input-start", "int"),
    end: getVal("input-end", "int"),
    ref: getVal("input-ref-base", "string", "G").toUpperCase(),
    alt: getVal("input-alt-base", "string", "A").toUpperCase(),
    strand: getVal("input-strand", "string", "+"),
    top_k: 5,

    // Amplicon
    min_amplicon_length: getVal("min_amplicon_length", "int", 60),
    max_amplicon_length: getVal("max_amplicon_length", "int", 150),

    // Primer
    primer_min_length: getVal("primer_min_length", "int", 20),
    primer_opt_length: getVal("primer_opt_length", "int", 25),
    primer_max_length: getVal("primer_max_length", "int", 30),

    primer_min_tm: getVal("primer_min_tm", "float", 55.0),
    primer_opt_tm: getVal("primer_opt_tm", "float", 60.0),
    primer_max_tm: getVal("primer_max_tm", "float", 65.0),

    primer_min_gc: getVal("primer_min_gc", "float", 35.0),
    primer_opt_gc: getVal("primer_opt_gc", "float", 50.0),
    primer_max_gc: getVal("primer_max_gc", "float", 65.0),

    // Probe
    probe_min_length: getVal("probe_min_length", "int", 20),
    probe_opt_length: getVal("probe_opt_length", "int", 25),
    probe_max_length: getVal("probe_max_length", "int", 30),

    min_primer_probe_tm_diff: getVal("min_primer_probe_tm_diff", "float", 5.0),
    max_primer_probe_tm_diff: getVal("max_primer_probe_tm_diff", "float", 10.0),

    // Absolute probe Tm
    probe_min_tm: getVal("probe_min_tm", "float", 65.0),
    probe_opt_tm: getVal("probe_opt_tm", "float", 67.0),
    probe_max_tm: getVal("probe_max_tm", "float", 70.0),

    probe_min_gc: getVal("probe_min_gc", "float", 35.0),
    probe_opt_gc: getVal("probe_opt_gc", "float", 50.0),
    probe_max_gc: getVal("probe_max_gc", "float", 65.0),

    // Constraints
    probe_max_poly_g: getVal("probe_max_poly_g", "int", 3),
    probe_max_3_end_gc: getVal("probe_max_3_end_gc", "int", 2),
    probe_avoid_5_prime_g: getVal("probe_avoid_5_prime_g", "bool", true),

    // QC
    qc_hairpin_min_dg: getVal("qc_hairpin_min_dg", "float", -6.0),
    qc_homodimer_min_dg: getVal("qc_homodimer_min_dg", "float", -6.0),
    qc_heterodimer_min_dg: getVal("qc_heterodimer_min_dg", "float", -6.0),
    qc_min_identity: getVal("qc_min_identity", "float", 90.0),
    qc_min_hit_length: getVal("qc_min_hit_length", "int", 13),
    qc_blast_max_alignments: getVal("qc_blast_max_alignments", "int", 50),
    qc_use_ispcr_check: getVal("qc_use_ispcr_check", "bool", false),
    qc_primer_max_diff_tm: getVal("qc_primer_max_diff_tm", "float", 3.0),

    // QC Probe (Constraints 재사용)
    qc_probe_avoid_5_prime_g: getVal("probe_avoid_5_prime_g", "bool", true),
    qc_probe_max_poly_g: getVal("probe_max_poly_g", "int", 3),
    qc_probe_max_3_end_gc: getVal("probe_max_3_end_gc", "int", 2),
  };
}

function validatePayload(p) {
  if (!p.chrom || !Number.isFinite(p.start)) {
    return "Please fill in Chromosome and Position.";
  }
  return null;
}

/* -----------------------------
   4) Render: Results / Metadata / Summary
------------------------------ */

function renderResults(results) {
  const tbody = byId("candidate-tbody");
  if (!tbody) return;

  tbody.innerHTML = "";

  const listData = results?.single_total_amplicons || [];
  if (!listData.length) {
    tbody.innerHTML =
      '<tr><td colspan="13" style="text-align:center; padding:40px;">No candidates found.</td></tr>';
    return;
  }

  listData.forEach((item, idx) => {
    const tr = document.createElement("tr");
    tr.style.cursor = "pointer";
    tr.id = `rank-row-${idx}`;

    tr.innerHTML = `
      <td class="text-center">${item.rank ?? (idx + 1)}</td>
      <td class="mono-cell">${item.forward_primer || "-"}</td>
      <td class="mono-cell">${item.reverse_primer || "-"}</td>
      <td class="mono-cell">${item.probe || "-"}</td>
      <td class="text-center">${formatNum1(item.tm_f)}</td>
      <td class="text-center">${formatNum1(item.tm_r)}</td>
      <td class="text-center">${formatNum1(item.tm_p)}</td>
      <td class="text-center">${formatNum1(item.gc_f)}%</td>
      <td class="text-center">${formatNum1(item.gc_r)}%</td>
      <td class="text-center">${formatNum1(item.gc_p)}%</td>
      <td class="text-center">${item.amplicon_size ?? 0}</td>
      <td class="text-nowrap">${item.genomic_pos || "-"}</td>
    `;

    tr.addEventListener("click", () => {
      document.querySelectorAll("#candidate-tbody tr").forEach((r) => r.classList.remove("selected-row"));
      tr.classList.add("selected-row");
      const alnView = byId("alignment-view");
      if (alnView) alnView.innerText = item.alignment_text_block || "No alignment data.";
    });

    tbody.appendChild(tr);
  });

  // 첫 행 자동 선택
  byId("rank-row-0")?.click();
}

function renderMetadata(results, payload) {
  const panel = byId("metadata-panel");
  const container = byId("metadata-table");
  if (!panel || !container || !results || !payload) return;

  panel.style.display = "block";

  const total = results.single_result?.total_count || 0;
  const filtered = results.single_total_amplicons?.length || 0;
  const rejected = total - filtered;

  container.innerHTML = `
    <div class="meta-report-grid">
      <div class="meta-card">
        <div class="card-header">📍 TARGET & AMPLICON</div>
        <div class="card-item"><span>Reference</span> <b>${payload.reference}</b></div>
        <div class="card-item"><span>Mutation</span> <b class="text-danger">${payload.ref} > ${payload.alt}</b></div>
        <div class="card-item"><span>Amp Size</span> <b>${payload.min_amplicon_length}-${payload.max_amplicon_length} bp</b></div>
        <div class="card-item"><span>Filtered</span> <b>${filtered}</b> / <b>${total}</b> <span style="opacity:.7">(rejected ${rejected})</span></div>
        <div class="help-text">
          * <b>Mutation:</b> Probe 서열 내부에 타겟 변이 부위를 포함하는 디자인인지 확인합니다.<br>
          * <b>Amp Size:</b> qPCR 서열(80-200bp)이 이상적입니다.
        </div>
      </div>

      <div class="meta-card">
        <div class="card-header">🔍 PROBE CONSTRAINTS</div>
        <div class="card-item"><span>Avoid 5' G</span> <b>${payload.qc_probe_avoid_5_prime_g ? "YES" : "NO"}</b></div>
        <div class="card-item"><span>Max Poly-G</span> <b>${payload.qc_probe_max_poly_g} nt</b></div>
        <div class="card-item"><span>3' End GC</span> <b>Max ${payload.qc_probe_max_3_end_gc}</b></div>
        <div class="help-text">
          * <b>Avoid 5' G:</b> 5' 말단의 G는 형광을 소광(Quenching)시켜 감도를 떨어뜨립니다.<br>
          * <b>Avoid 3' G-C Count &lt; 3 :</b> 결합력이 강한 GC 수가 적은 probe 선정으로 Non-specific binding 감소.<br>
          * <b>Poly-G:</b> 4개 이상의 연속된 G는 G-quadruplex 구조를 형성할 위험이 있습니다.
        </div>
      </div>

      <div class="meta-card">
        <div class="card-header">🛡️ BLAST & SPECIFICITY</div>
        <div class="card-item"><span>Min Identity</span> <b>${payload.qc_min_identity}%</b></div>
        <div class="card-item"><span>Min Hit BP</span> <b>${payload.qc_min_hit_length} bp</b></div>
        <div class="card-item"><span>Max Align</span> <b>${payload.qc_blast_max_alignments} Hits</b></div>
        <div class="help-text">
          * <b>Min Identity:</b> 유전체 내 다른 부위와의 상동성 허용치입니다. 높을수록 엄격합니다.<br>
          * <b>Min Hit BP:</b> 연속적으로 일치한 염기서열 수.<br>
          * <b>Max Align:</b> Blast 결과 중 상위 수의 위치만 비특이적 결합 유무 확인.
        </div>
      </div>
    </div>
  `;
}

function renderSummary(payload, results) {
  const now = new Date();
  const dateStr = now.toLocaleDateString() + " " + now.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });

  setText("summary-chrom", payload.chrom);
  setText("summary-pos", `${payload.start} - ${payload.end}`);
  setText("summary-date", dateStr);
  setText("summary-reference", payload.reference);
  setText("summary-mutation", `${payload.ref} > ${payload.alt}`);

  if (results?.status === "success") setStatus("COMPLETED", "green");
  else if (results?.status) setStatus("FAILED", "red");
  else setStatus("READY", "#2b384c");
}

/* -----------------------------
   5) Export: JSON / HTML
------------------------------ */

function exportJSON() {
  if (!state.lastResults) {
    alert("내보낼 데이터가 없습니다. 먼저 분석을 실행하세요.");
    return;
  }

  const projectOutput = {
    project_info: {
      name: `qPCR_Project_${state.lastResults.export_meta?.region?.name || "Design"}`,
      date: new Date().toISOString(),
      app_version: "1.0.0",
    },
    input_parameters: state.lastPayload || {},
    design_results: state.lastResults,
  };

  const jsonString = JSON.stringify(projectOutput, null, 4);
  const blob = new Blob([jsonString], { type: "application/json" });
  const url = URL.createObjectURL(blob);

  const link = document.createElement("a");
  link.href = url;
  link.download = `qPCR_RawData_${Date.now()}.json`;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}

function exportHTML() {
  // main export 영역: .main-panel 우선, 없으면 .main
  const target = $(".main-panel") || $(".main");
  if (!target) {
    alert("리포트 영역(.main-panel/.main)을 찾을 수 없습니다.");
    return;
  }

  // 외부 도메인 CSSRules 접근은 실패할 수 있음(try/catch)
  const styles = Array.from(document.styleSheets)
    .map((ss) => {
      try {
        return Array.from(ss.cssRules).map((r) => r.cssText).join("\n");
      } catch {
        return "";
      }
    })
    .join("\n");

  const reportArea = target.innerHTML;

  const htmlContent = `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <title>qPCR Design Report - ${new Date().toLocaleString()}</title>
  <style>
    body { font-family: 'Nunito', sans-serif; padding: 20px; background: #f5f5f5; }
    ${styles}
    .no-print, .export-buttons { display: none !important; }
    .report-wrap { background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
    .table-container { overflow-x: auto !important; }
  </style>
</head>
<body>
  <div class="report-wrap">
    ${reportArea}
  </div>
</body>
</html>`;

  const blob = new Blob([htmlContent], { type: "text/html;charset=utf-8;" });
  const url = URL.createObjectURL(blob);

  const link = document.createElement("a");
  link.href = url;
  link.download = `qPCR_Report_${Date.now()}.html`;
  link.click();
  URL.revokeObjectURL(url);
}

/* -----------------------------
   6) Main Run (Fetch)
------------------------------ */

async function runDesign() {
  const runBtn = byId("btn-run-qpcr");
  if (!runBtn) return;

  // 저장 + payload
  saveInputsToStorage();
  const payload = buildPayload();
  const err = validatePayload(payload);
  if (err) {
    alert(err);
    return;
  }

  // UI: loading
  const originalText = runBtn.innerText;
  runBtn.innerText = "⏳ RUNNING...";
  runBtn.disabled = true;
  setStatus("RUNNING", "#2b384c");

  try {
    console.log("📤 Sending payload:", payload);

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
    console.log("📥 Received results:", results);

    // state 저장
    state.lastPayload = payload;
    state.lastResults = results;

    // 렌더
    renderResults(results);
    renderMetadata(results, payload);
    renderSummary(payload, results);

  } catch (e) {
    console.error("Fetch Error:", e);
    alert("Design failed: " + e.message);
    setStatus("ERROR", "red");
  } finally {
    runBtn.innerText = originalText;
    runBtn.disabled = false;
  }
}

/* -----------------------------
   7) Wire up Events (single entry)
------------------------------ */

document.addEventListener("DOMContentLoaded", () => {
  // restore inputs
  restoreInputsFromStorage();

  // run
  byId("btn-run-qpcr")?.addEventListener("click", (e) => {
    e.preventDefault();
    runDesign();
  });

  // exports
  byId("btn-export-json")?.addEventListener("click", exportJSON);
  byId("btn-export-html")?.addEventListener("click", exportHTML);

  console.log("✅ qPCR app ready.");
});