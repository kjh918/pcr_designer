/**
 * app.js (AS-PCR Accordion UI & Unified Params Version)
 * - LocalStorage 자동 저장/복구
 * - FastAPI 통신 (단일 서열 'sequence' 파라미터 기반)
 * - AS-PCR Set 구조(WT, ALT, WT_MM, ALT_MM) 완벽 병합(Accordion) 렌더링
 */

const state = {
    lastResults: null,
    lastPayload: null,
};

const STORAGE_PREFIX = "aspcr_";

// 🔥 공통 params.html에 있는 실제 ID(mspcr-*)들을 완벽하게 매핑합니다.
const TRACKED_INPUTS = [
    "mspcr-design-name", "aspcr-design-name", "input-design-name", 
    "mspcr-raw-sequence", "aspcr-raw-sequence", "manual-sequence", 
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

const API_URL = "http://192.168.0.35:9000/api/design/aspcr";

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

function extractMutation(seq) {
    const match = seq.match(/\[(.*?)\]/);
    if (match) {
        const alleles = match[1].replace(/\s/g, '').replace('/', ',');
        const parts = alleles.split(',');
        if (parts.length >= 2) return `${parts[0]} > ${parts[1]}`;
        return alleles;
    }
    return "-";
}

function buildPayload() {
    // 🔥 공용 params.html의 실제 ID(mspcr-raw-sequence)를 최우선으로 찾습니다.
    const rawSeq = (byId("mspcr-raw-sequence")?.value || byId("aspcr-raw-sequence")?.value || byId("manual-sequence")?.value || "").replace(/\s/g, "").toUpperCase();
    const designName = getVal("mspcr-design-name", "string", getVal("aspcr-design-name", "string", getVal("input-design-name", "string", "ASPCR_Design")));

    return {
        design_name: designName,
        sequence: rawSeq,  
        reference_genome: getVal("input-ref-genome", "string", "hg38"),
        top_k: 5,
        
        fixed_prime: getVal("aspcr_fixed_prime", "string", "forward"),
        mismatch_pos: getVal("aspcr_mismatch_pos", "int", 3),
        mismatch_intensity: getVal("aspcr_mismatch_intensity", "string", "strong"),
        
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
                max_tm_diff: getVal("qc_primer_max_diff_tm", "float", 5.0)
            }
        }
    };
}

// 🔥 VCF 좌표 검증을 제거하고, 대괄호 포함 여부만 검증합니다.
function validatePayload(p) {
    if (!p.sequence || p.sequence.length < 30) {
        return "Please enter a valid DNA sequence (at least 30bp).";
    }
    if (!p.sequence.includes("[") || !p.sequence.includes("]")) {
        return "서열에 대괄호('[', ']')를 사용하여 변이 타겟을 명시해야 합니다. (e.g. ATGC[A,G]ATGC)";
    }
    return null;
}

function renderSummary(payload, results) {
    const summary = results?.summary || {};
    const now = new Date();
    const dateStr = now.toLocaleDateString() + " " + now.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });

    setText("summary-chrom", "Auto");
    setText("summary-pos", "Auto");
    setText("summary-date", dateStr);
    setText("summary-reference", payload.reference_genome || "hg38");
    setText("summary-mutation", extractMutation(payload.sequence));

    if (results?.status === "success") {
        const passedCount = summary.passed_count || 0;
        const totalCount = summary.total_count || 0;
        const allFailed = passedCount === 0 && totalCount > 0;
        setStatus(allFailed ? `FAIL (${passedCount} PASS)` : "SUCCESS", "color");
    } else if (results?.status) {
        setStatus("ERROR", "color");
    } else {
        setStatus("READY", "#2b384c");
    }
}

// 🔥 [아코디언 토글 전역 함수]
window.toggleSet = function(setId) {
    const rows = document.querySelectorAll('.' + setId);
    const icon = document.getElementById('icon-' + setId);
    let isHidden = true;
    
    rows.forEach(row => {
        if (row.style.display === 'none') {
            row.style.display = 'table-row';
            isHidden = false;
        } else {
            row.style.display = 'none';
        }
    });
    
    if (icon) {
        icon.textContent = isHidden ? '+' : '-';
    }
};

// 🔥 [아코디언 렌더러] 17열 구조에 맞춰 메인(Set) 행 클릭 시 서브(Allele) 행들이 토글됩니다.
function renderResults(results) {
    const tbody = document.getElementById("candidate-tbody");
    const alignmentView = document.getElementById("alignment-view");
    
    if (!tbody) return;

    const ampliconsData = results?.results || [];
    
    if (!ampliconsData || ampliconsData.length === 0) {
        tbody.innerHTML = `<tr><td colspan="17" style="text-align:center; padding:60px; color:#999; font-size: 13px;">No candidates found or all failed.</td></tr>`;
        return;
    }

    let html = "";
    let alignmentTexts = [];

    ampliconsData.forEach((set, setIndex) => {
        const rank = set.rank || (setIndex + 1);
        const setId = set.set_id || `Set_${rank}`;
        const setPass = set.set_qc_pass;
        
        const sampleAllele = Object.values(set.alleles || {})[0];
        if (sampleAllele && sampleAllele.amplicon_info?.alignment_text_block) {
            alignmentTexts.push(`[${setId} ALIGNMENT]\n${sampleAllele.amplicon_info.alignment_text_block}`);
        }

        const qcBadge = setPass ? `<span style="color:#198754; font-weight:bold; font-size:12px;">PASS</span>` : `<span style="color:#d9534f; font-weight:bold; font-size:12px;">FAIL</span>`;
        const setRowId = `set-row-${rank}`;
        
        // 1. 요약 메인 행 (Set 정보) - 클릭 시 toggleSet 호출
        html += `
            <tr style="background-color: #fdfdfe; cursor: pointer; border-top: 2px solid #adb5bd;" onclick="toggleSet('${setRowId}')">
                <td style="text-align: center; border-bottom: 1px solid #dee2e6;">
                    <span id="icon-${setRowId}" style="font-size: 16px; font-weight: bold; color: #6c757d; border: 1px solid #ccc; padding: 0 5px; border-radius: 3px;">+</span>
                </td>
                <td style="font-weight: bold; border-bottom: 1px solid #dee2e6; color: #495057;">Rank ${rank}</td>
                <td style="text-align: center; font-weight: bold; border-bottom: 1px solid #dee2e6; border-right: 2px solid #adb5bd; color: #333;">${setId}<br>${qcBadge}</td>
                <td colspan="14" style="border-bottom: 1px solid #dee2e6; color: #6c757d; font-size: 12px; vertical-align: middle;">
                    👉 <em>Click here to expand and view details for <strong>WT, ALT, WT_MM, ALT_MM</strong> alleles.</em>
                </td>
            </tr>
        `;

        // 2. 상세 서브 행 (각 Allele의 Fwd, Rev, HD, QC 정보)
        const alleleKeys = ["wt", "alt", "wt_mm", "alt_mm"];
        alleleKeys.forEach((key, index) => {
            const alleleData = set.alleles[key];
            if (!alleleData) return;

            const fwd = alleleData.oligos?.forward || { sequence: "-", tm: 0, gc: 0 };
            const rev = alleleData.oligos?.reverse || { sequence: "-", tm: 0, gc: 0 };
            const hetero = alleleData.oligos?.heterodimer || { fr_dg: 0, fp_dg: 0, rp_dg: 0 };
            
            const isIndivPass = alleleData.qc_info?.is_pass;
            const qcText = isIndivPass ? "PASS" : (alleleData.qc_info?.fail_reason || "FAIL");
            const qcColor = isIndivPass ? "#198754" : "#d9534f";
            
            let labelColor = key.includes("mm") ? "#d9534f" : "#0275d8"; 
            let rowBottomBorder = (index === 3) ? `border-bottom: 2px solid #adb5bd;` : `border-bottom: 1px solid #e9ecef;`;

            html += `
            <tr class="allele-detail-row ${setRowId}" style="display: none; background-color: #ffffff;">
                <td style="${rowBottomBorder}"></td>
                <td style="font-weight: bold; color: ${labelColor}; ${rowBottomBorder}">${key.toUpperCase()}</td>
                <td style="${rowBottomBorder} border-right: 2px solid #adb5bd;"></td>
                
                <!-- Fwd (5 cols) -->
                <td class="mono-cell" style="${rowBottomBorder}">${fwd.sequence}</td>
                <td class="text-center" style="${rowBottomBorder}">${formatNum1(fwd.tm)}</td>
                <td class="text-center" style="${rowBottomBorder}">${formatNum1(fwd.gc)}</td>
                <td class="text-center" style="${rowBottomBorder} color: #d9534f;">${formatNum1(fwd.hairpin_dg)}</td>
                <td class="text-center" style="${rowBottomBorder} color: #d9534f; border-right: 2px solid #adb5bd;">${formatNum1(fwd.homodimer_dg)}</td>

                <!-- Rev (5 cols) -->
                <td class="mono-cell" style="${rowBottomBorder}">${rev.sequence}</td>
                <td class="text-center" style="${rowBottomBorder}">${formatNum1(rev.tm)}</td>
                <td class="text-center" style="${rowBottomBorder}">${formatNum1(rev.gc)}</td>
                <td class="text-center" style="${rowBottomBorder} color: #d9534f;">${formatNum1(rev.hairpin_dg)}</td>
                <td class="text-center" style="${rowBottomBorder} color: #d9534f; border-right: 2px solid #adb5bd;">${formatNum1(rev.homodimer_dg)}</td>

                <!-- Hetero (1 col) -->
                <td class="text-center" style="${rowBottomBorder} color: #d9534f; border-right: 2px solid #adb5bd;">${formatNum1(hetero.fr_dg)}</td>

                <!-- Amp/QC (3 cols) -->
                <td class="text-center" style="${rowBottomBorder}">${formatNum1(alleleData.metrics?.pair_penalty)}</td>
                <td class="text-center text-nowrap" style="${rowBottomBorder}">${alleleData.amplicon_info?.genomic_pos || "-"}</td>
                <td style="${rowBottomBorder} color: ${qcColor}; font-size: 10px; max-width: 200px; word-wrap: break-word;">${qcText}</td>
            </tr>
            `;
        });
    });

    tbody.innerHTML = html;

    if (alignmentView) {
        if (alignmentTexts.length > 0) {
            alignmentView.textContent = alignmentTexts.join("\n\n" + "=".repeat(100) + "\n\n");
        } else {
            alignmentView.textContent = "No alignment data available.";
        }
    }
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
    const projectName = projectOutput.metadata?.project_name || "ASPCR_Project";

    const jsonString = JSON.stringify(projectOutput, null, 4);
    const blob = new Blob([jsonString], { type: "application/json" });
    const url = URL.createObjectURL(blob);

    const link = document.createElement("a");
    link.href = url;
    link.download = `${projectName}_Report.json`;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
}

function exportHTML(e) {
    if (e) e.preventDefault();
    if (!state.lastResults) {
        alert("내보낼 데이터가 없습니다.");
        return;
    }
    const target = document.querySelector(".main-content");
    if (!target) return;

    let styles = "";
    try {
        styles = Array.from(document.styleSheets)
            .map((ss) => {
                try { return Array.from(ss.cssRules).map((r) => r.cssText).join("\n"); } 
                catch { return ""; }
            }).join("\n");
    } catch(err) {}

    const projectName = state.lastPayload?.design_name || "ASPCR_Project";
    const htmlContent = `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <title>AS-PCR Design Report - ${projectName}</title>
  <style>
    body { font-family: 'Nunito', sans-serif; padding: 20px; background: #f5f5f5; }
    ${styles}
    .no-print, details { display: none !important; }
    .panel { background: white; padding: 20px; border-radius: 8px; border: 1px solid #ddd; margin-bottom: 20px; }
  </style>
</head>
<body>
  <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1);">
    <h2>🔍 AS-PCR Design Report</h2>
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
    const runBtn = byId("btn-run-aspcr");
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
        console.log("📤 Sending AS-PCR Payload:", payload);
        const res = await fetch(API_URL, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(payload),
        });

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
        console.log("📥 Received Results:", results);

        state.lastPayload = payload;
        state.lastResults = results;

        if (results.status !== "success" || results.error || results.reason) {
            const errContent = byId("error-log-content");
            if (errPanel && errContent) {
                errPanel.style.display = "block";
                let logs = results.log_messages || results.error || results.reason || "Unknown error occurred.";
                errContent.innerText = Array.isArray(logs) ? logs.join('\n') : logs;
            }
        }

        renderResults(results);
        renderSummary(payload, results);

        // 첫 번째 행은 자동으로 펼쳐줍니다
        setTimeout(() => window.toggleSet("set-row-1"), 100);

    } catch (e) {
        console.error("Fetch Error:", e);
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

/* -----------------------------
   7) Wire up Events
------------------------------ */

document.addEventListener("DOMContentLoaded", () => {
    restoreInputsFromStorage();

    byId("btn-run-aspcr")?.addEventListener("click", (e) => {
        e.preventDefault();
        runDesign();
    });

    byId("btn-export-json")?.addEventListener("click", exportJSON);
    byId("btn-export-html")?.addEventListener("click", exportHTML);

    byId("btn-load-tutorial")?.addEventListener("click", (e) => {
        e.preventDefault();
        const designNameInput = byId("mspcr-design-name") || byId("aspcr-design-name");
        if (designNameInput) designNameInput.value = "EGFR_L858R_Tutorial";
        
        const seqInput = byId("mspcr-raw-sequence") || byId("aspcr-raw-sequence");
        if (seqInput) seqInput.value = "GAAAATGACAAAGAACAGCTCAAAGCAATTTCTACACGAGATCCTCTCTCTGAAATCACT[G,A]AGCAGGAGAAAGATTTTCTATGGAGTCACAGGTAAGTGCTA";
        
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