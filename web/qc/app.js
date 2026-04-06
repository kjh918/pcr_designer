/**
 * qc_app.js (Complete Version)
 * - Payload 평탄화(Flat) 적용 및 최적화된 JSON 응답(DRY 패턴) 파싱 완료
 * - 26-column Detailed Analytics Table 적용 완료
 * - Probe 및 Primer 분리된 QC Parameters 파싱 완료
 */

const state = {
    lastResults: null,
    lastPayload: null,
};

const STORAGE_PREFIX = "qc_";

const TRACKED_INPUTS = [
    "qc-project-name", "qc-forward-seq", "qc-reverse-seq", "qc-probe-seq", "qc-template-seq",
    "qc-ref-genome",
    "qc_hairpin_min_dg", "qc_homodimer_min_dg", "qc_heterodimer_min_dg",
    "qc_min_identity", "qc_min_hit_length", "qc_blast_max_alignments",
    "qc_min_amp_size", "qc_max_amp_size", "qc_use_ispcr_check", 
    "qc_primer_min_diff_tm", "qc_primer_max_diff_tm",
    "qc_probe_min_tm_diff", "qc_probe_max_tm_diff", "qc_probe_poly_g", "qc_probe_avoid_5g"
];

// 🔥 로컬 환경 포트 확인
const API_URL = "http://192.168.0.35:9000/api/design/qc";

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
    if (type === "int") return Number.isFinite(parseInt(raw, 10)) ? parseInt(raw, 10) : def;
    if (type === "float") return Number.isFinite(parseFloat(raw)) ? parseFloat(raw) : def;
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

const cleanSequence = (str) => str ? str.replace(/[\s0-9]/g, '').toUpperCase() : "";

function saveInputsToStorage() {
    TRACKED_INPUTS.forEach((id) => {
        const el = byId(id);
        if (!el) return;
        localStorage.setItem(STORAGE_PREFIX + id, String(el.type === "checkbox" ? el.checked : el.value));
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

function buildPayload() {
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
            hairpin_min_dg: getVal("qc_hairpin_min_dg", "float", -5.0),
            homodimer_min_dg: getVal("qc_homodimer_min_dg", "float", -6.0),
            heterodimer_min_dg: getVal("qc_heterodimer_min_dg", "float", -6.0),
            min_identity: getVal("qc_min_identity", "float", 90.0),
            min_hit_length: getVal("qc_min_hit_length", "int", 13),
            blast_max_alignments: getVal("qc_blast_max_alignments", "int", 50),
            min_amp_size: getVal("qc_min_amp_size", "int", 50),
            max_amp_size: getVal("qc_max_amp_size", "int", 300),
            use_ispcr_check: getVal("qc_use_ispcr_check", "bool", false),
            
            // 🔥 추가된 프라이머 및 프로브 조건 패키징
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
    if (!p.sequences.forward || !p.sequences.reverse) return "Forward Primer와 Reverse Primer는 필수 입력 항목입니다.";
    if (p.sequences.forward.length < 15 || p.sequences.reverse.length < 15) return "프라이머 길이는 최소 15bp 이상이어야 합니다.";
    return null;
}

// 🔥 핵심 서버 통신 함수
async function runDesign() {
    const runBtn = byId("btn-run-qc");
    if (!runBtn) return;

    saveInputsToStorage();
    const payload = buildPayload();
    const err = validatePayload(payload);
    if (err) { alert(err); return; }

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
        renderSummary(results);

    } catch (e) {
        alert("QC Request failed: " + e.message);
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

// ---------------------------------------------------------
// 렌더링(화면 출력) 함수들
// ---------------------------------------------------------
function renderSummary(results) {
    // 🔥 백엔드 JSON 최적화 구조 반영 (results.inputs, results.metadata)
    const meta = results.metadata;
    const inputs = results.inputs || {};
    const seqs = inputs.sequences || {};
    const qc = inputs.qc_criteria || {};
    const summary = results.summary || {};

    if (!meta) return;

    setText("summary-project", meta.project_name || "QC_Project");
    setText("summary-reference", meta.reference_genome || "None");
    
    const dateObj = meta.timestamp ? new Date(meta.timestamp) : new Date();
    setText("summary-date", dateObj.toLocaleDateString() + " " + dateObj.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" }));

    const passedCount = summary.passed_count || 0;
    const totalCount = summary.total_count || 0;
    setText("summary-counts", `${passedCount} / ${totalCount}`);

    if (results.status === "success") {
        const allFailed = passedCount === 0 && totalCount > 0;
        setStatus(allFailed ? `FAIL (${passedCount} PASS)` : "SUCCESS", "color");
    } else {
        setStatus("ERROR", "color");
    }

    setText("summary-seq-fwd", seqs.forward || "-");
    setText("summary-seq-rev", seqs.reverse || "-");
    setText("summary-seq-prb", seqs.probe || "N/A");

    // Thermodynamics
    const thermo = qc.thermodynamics || qc;
    setText("summary-qc-hairpin", `${thermo.hairpin_min_dg ?? "-5"} kcal`);
    setText("summary-qc-homo", `${thermo.homodimer_min_dg ?? "-6"} kcal`);
    setText("summary-qc-hetero", `${thermo.heterodimer_min_dg ?? "-6"} kcal`);
    
    // Primer QC
    const minTmDiff = qc.primer?.min_diff_tm ?? "0";
    const maxTmDiff = qc.primer?.max_diff_tm ?? qc.max_tm_diff ?? "3";
    setText("summary-qc-tmdiff", `${minTmDiff}~${maxTmDiff} ℃`);

    // Probe QC
    const probeQc = qc.probe || {};
    setText("summary-qc-probe-5g", probeQc.avoid_5_prime_g === false ? "Allow" : "Avoid");
    setText("summary-qc-probe-polyg", `≤ ${probeQc.max_probe_poly_g ?? "3"}`);
    setText("summary-qc-probe-tmdiff", `${probeQc.min_primer_probe_tm_diff ?? "5"}~${probeQc.max_primer_probe_tm_diff ?? "10"} ℃`);

    // Specificity & Size
    const blast = qc.blast || qc;
    setText("summary-qc-ident", `${blast.min_identity ?? "90"} %`);
    setText("summary-qc-hitlen", `${blast.min_hit_length ?? "13"} bp`);
    setText("summary-qc-aligns", `${blast.blast_max_alignments ?? blast.max_alignments ?? "50"}`);

    const amp = qc.amplicon || qc;
    const minS = amp.min_amp_size ?? amp.min_size ?? "50";
    const maxS = amp.max_amp_size ?? amp.max_size ?? "300";
    setText("summary-qc-ampsize", `${minS}~${maxS} bp`);
}

function renderResults(results) {
    const tbody = document.getElementById("candidate-tbody");
    if (!tbody) return;
    tbody.innerHTML = "";

    const listData = results?.results || [];
    if (!listData.length) {
        // 🔥 26컬럼 콜스팬 적용
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

        // 🔥 백엔드 Schema에서 던져주는 객체 완벽 매핑
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

        // 🔥 정확히 26개의 <td> 태그 생성
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
            const alnView = document.getElementById("alignment-view");
            if (alnView) alnView.innerText = alignText;
        });

        tbody.appendChild(tr);
    });

    // 자동 첫 행 클릭
    setTimeout(() => {
        const firstRow = document.getElementById("rank-row-0");
        if (firstRow) firstRow.click();
    }, 50);
}

// ---------------------------------------------------------
// Export 함수들
// ---------------------------------------------------------
function exportJSON(e) {
    if (e) e.preventDefault();
    if (!state.lastResults) {
        alert("내보낼 데이터가 없습니다. 먼저 QC 분석을 실행하세요.");
        return;
    }

    // 🔥 프론트엔드에서 더 이상 억지로 감싸지 않고, 백엔드에서 만든 완벽한 구조를 그대로 사용합니다.
    const projectOutput = state.lastResults;
    const projectName = projectOutput.metadata?.project_name || "QC_Project";

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
        alert("내보낼 데이터가 없습니다. 먼저 QC 분석을 실행하세요.");
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

    const projectName = state.lastPayload?.project_name || "QC_Project";
    const htmlContent = `<!DOCTYPE html>
<html lang="ko">
<head>
  <meta charset="UTF-8" />
  <title>QC Report - ${projectName}</title>
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
    <h2 style="color: #1e707a; margin-top: 0;">🔍 Sequence QC Analysis Report</h2>
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
    link.download = `${projectName}_QC_Report.html`;
    link.click();
    URL.revokeObjectURL(url);
}

// ---------------------------------------------------------
// 이벤트 리스너 세팅
// ---------------------------------------------------------
document.addEventListener("DOMContentLoaded", () => {
    restoreInputsFromStorage();

    byId("btn-run-qc")?.addEventListener("click", (e) => { e.preventDefault(); runDesign(); });
    
    byId("btn-export-json")?.addEventListener("click", exportJSON);
    byId("btn-export-html")?.addEventListener("click", exportHTML);

    byId("btn-load-tutorial")?.addEventListener("click", (e) => {
        e.preventDefault();
        const pName = byId("qc-project-name");
        const fwd = byId("qc-forward-seq");
        const rev = byId("qc-reverse-seq");
        const prb = byId("qc-probe-seq");
        const ref = byId("qc-ref-genome");
        
        if (pName) pName.value = "EGFR_L858R_QC_Test";
        if (fwd) fwd.value = "GCACATCGAGGCCAACACT";
        if (rev) rev.value = "TGAGAAAAGTGAACCTGCGGA";
        if (prb) prb.value = "TTGCCCGGACCTGGGATGC";
        if (ref) ref.value = "hg38";
        
        const btn = e.target;
        btn.innerText = "✅ Loaded!";
        btn.style.backgroundColor = "#28a745";
        setTimeout(() => {
            btn.innerText = "🧪 Tutorial";
            btn.style.backgroundColor = "#17a2b8";
        }, 2000);
    });

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

    document.querySelectorAll(".seq-input").forEach(input => {
        input.addEventListener("paste", function(e) {
            e.preventDefault();
            const text = (e.clipboardData || window.clipboardData).getData("text");
            this.value = cleanSequence(text);
        });
    });

    byId("btn-clear-seqs")?.addEventListener("click", function(e) {
        e.preventDefault();
        document.querySelectorAll(".seq-input").forEach(input => input.value = "");
        if (byId("qc-project-name")) byId("qc-project-name").value = "";
    });
});