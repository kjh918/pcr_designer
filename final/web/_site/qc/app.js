/**
 * qc_app.js (Refactored)
 * - LocalStorage 자동 저장/복구
 * - Quick Paste (엑셀 복붙) 및 서열 정제
 * - FastAPI 통신 (QC 전용 Payload)
 * - 결과/메타 렌더링 및 Export (JSON / HTML)
 */

/* -----------------------------
   0) App State / Constants
------------------------------ */
const state = {
    lastResults: null,
    lastPayload: null,
};

const STORAGE_PREFIX = "qc_";

// 🔥 HTML 폼의 ID와 정확히 일치하도록 수정
const TRACKED_INPUTS = [
    "qc-project-name", "qc-forward-seq", "qc-reverse-seq", "qc-probe-seq", "qc-template-seq",
    "qc-ref-genome",
    "qc_hairpin_min_dg", "qc_homodimer_min_dg", "qc_heterodimer_min_dg",
    "qc_min_identity", "qc_min_hit_length", "qc_blast_max_alignments",
    "qc_min_amp_size", "qc_max_amp_size", "qc_use_ispcr_check", 
    "qc_primer_max_diff_tm"
];

const API_URL = "http://192.168.0.35:9000/api/design/qc";

/* -----------------------------
   1) Utils (DOM / format / read)
------------------------------ */
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
        return type === "bool" ? !!el.checked : el.checked;
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

// 공백, 줄바꿈, 숫자 자동 제거 필터
const cleanSequence = (str) => {
    return str ? str.replace(/[\s0-9]/g, '').toUpperCase() : "";
};

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
    // 🔥 QC 전용 JSON 구조로 매핑
    return {
        project_name: getVal("qc-project-name", "string", "QC_Project"),
        sequences: {
            forward: cleanSequence(getVal("qc-forward-seq")),
            reverse: cleanSequence(getVal("qc-reverse-seq")),
            probe: cleanSequence(getVal("qc-probe-seq")),
            template: cleanSequence(getVal("qc-template-seq"))
        },
        reference_genome: getVal("qc-ref-genome", "string", "none"),
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
    if (!p.sequences.forward || !p.sequences.reverse) {
        return "Forward Primer와 Reverse Primer는 필수 입력 항목입니다.";
    }
    if (p.sequences.forward.length < 15 || p.sequences.reverse.length < 15) {
        return "프라이머 길이는 최소 15bp 이상이어야 합니다.";
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
        tbody.innerHTML = '<tr><td colspan="13" style="text-align:center; padding:40px;">No candidates found.</td></tr>';
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

    byId("rank-row-0")?.click();
}

function renderMetadata(results, payload) {
    const panel = byId("metadata-panel");
    const container = byId("metadata-table");
    if (!panel || !container || !results || !payload) return;

    panel.style.display = "block";

    container.innerHTML = `
        <div class="meta-report-grid">
            <div class="meta-card">
                <div class="card-header">📍 INPUT SEQUENCES</div>
                <div class="card-item"><span>Fwd Primer</span> <b>${payload.sequences.forward}</b></div>
                <div class="card-item"><span>Rev Primer</span> <b>${payload.sequences.reverse}</b></div>
                <div class="card-item"><span>Probe</span> <b>${payload.sequences.probe || "N/A"}</b></div>
            </div>
            <div class="meta-card">
                <div class="card-header">🛡️ THERMODYNAMICS & BLAST</div>
                <div class="card-item"><span>Hairpin dG</span> <b>Max ${payload.qc_criteria.thermodynamics.hairpin_min_dg}</b></div>
                <div class="card-item"><span>Dimer dG</span> <b>Max ${payload.qc_criteria.thermodynamics.homodimer_min_dg}</b></div>
                <div class="card-item"><span>BLAST Hit</span> <b>Min ${payload.qc_criteria.blast.min_hit_length}bp / ${payload.qc_criteria.blast.min_identity}%</b></div>
            </div>
        </div>
    `;
}

function renderSummary(payload, results) {
    const now = new Date();
    const dateStr = now.toLocaleDateString() + " " + now.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });

    setText("summary-date", dateStr);
    setText("summary-reference", payload.reference_genome);

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
            name: state.lastPayload?.project_name || "QC_Project",
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
    link.download = `${projectOutput.project_info.name}_RawData.json`;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
}

function exportHTML() {
    const target = $(".main-content") || $(".main");
    if (!target) {
        alert("리포트 영역을 찾을 수 없습니다.");
        return;
    }

    const styles = Array.from(document.styleSheets)
        .map((ss) => {
            try { return Array.from(ss.cssRules).map((r) => r.cssText).join("\n"); } 
            catch { return ""; }
        }).join("\n");

    const reportArea = target.innerHTML;

    const htmlContent = `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <title>QC Report - ${new Date().toLocaleString()}</title>
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
    link.download = `QC_Report_${Date.now()}.html`;
    link.click();
    URL.revokeObjectURL(url);
}

/* -----------------------------
   6) Main Run (Fetch)
------------------------------ */
async function runDesign() {
    const runBtn = byId("btn-run-qc");
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

        state.lastPayload = payload;
        state.lastResults = results;

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
   7) Init & Event Listeners
------------------------------ */
document.addEventListener("DOMContentLoaded", () => {
    // 1. Storage 복구
    restoreInputsFromStorage();

    // 2. 버튼 이벤트 연결
    byId("btn-run-qc")?.addEventListener("click", (e) => {
        e.preventDefault();
        runDesign();
    });
    byId("btn-export-json")?.addEventListener("click", exportJSON);
    byId("btn-export-html")?.addEventListener("click", exportHTML);

    // 3. Quick Paste (엑셀 복붙 처리기)
    const quickPasteBox = byId("qc-quick-paste");
    if (quickPasteBox) {
        quickPasteBox.addEventListener("paste", function(e) {
            e.preventDefault();
            const text = (e.clipboardData || window.clipboardData).getData("text");
            const parts = text.split(/[\t\n]+/).map(s => s.trim()).filter(s => s.length > 0);
            
            if (parts[0]) byId("qc-forward-seq").value = cleanSequence(parts[0]);
            if (parts[1]) byId("qc-reverse-seq").value = cleanSequence(parts[1]);
            if (parts[2]) byId("qc-probe-seq").value = cleanSequence(parts[2]);
            if (parts[3]) byId("qc-template-seq").value = cleanSequence(parts[3]);

            quickPasteBox.value = "✅ 서열이 성공적으로 분류 및 입력되었습니다!";
            quickPasteBox.style.color = "#28a745";
            quickPasteBox.style.fontWeight = "bold";
            
            setTimeout(() => {
                quickPasteBox.value = "";
                quickPasteBox.style.color = "";
                quickPasteBox.style.fontWeight = "normal";
            }, 2000);
        });
    }

    // 4. 개별 서열 입력창 자동 정제 (붙여넣기 시)
    document.querySelectorAll(".seq-input").forEach(input => {
        input.addEventListener("paste", function(e) {
            e.preventDefault();
            const text = (e.clipboardData || window.clipboardData).getData("text");
            this.value = cleanSequence(text);
        });
    });

    // 5. 서열 초기화 버튼
    byId("btn-clear-seqs")?.addEventListener("click", function(e) {
        e.preventDefault();
        document.querySelectorAll(".seq-input").forEach(input => input.value = "");
        byId("qc-project-name").value = "";
    });

    console.log("✅ QC app initialized successfully.");
});