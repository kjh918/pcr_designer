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

function renderResults(results) {
    const tbody = document.getElementById("candidate-tbody");
    if (!tbody) return;
    
    const ampliconsData = results?.results || [];
    if (!ampliconsData.length) {
        tbody.innerHTML = `<tr><td colspan="17" style="text-align:center; padding:60px; color:#999;">No candidates found.</td></tr>`;
        return;
    }

    let html = "";
    ampliconsData.forEach((set, idx) => {
        const rank = set.rank || (idx + 1);
        const setId = set.set_id || `Set_${rank}`;
        const setRowId = `set-row-${rank}`;
        const qcBadge = set.set_qc_pass ? `<span style="color:#198754; font-weight:bold;">PASS</span>` : `<span style="color:#d9534f; font-weight:bold;">FAIL</span>`;

        html += `<tr style="background-color: #fdfdfe; cursor: pointer; border-top: 2px solid #adb5bd;" onclick="toggleSet('${setRowId}')">
                <td style="text-align: center;"><span id="icon-${setRowId}" style="font-size: 16px; border: 1px solid #ccc; padding: 0 5px; border-radius: 3px;">+</span></td>
                <td style="font-weight: bold;">Rank ${rank}</td>
                <td style="text-align: center; font-weight: bold; border-right: 2px solid #adb5bd;">${setId}<br>${qcBadge}</td>
                <td colspan="14" style="color: #6c757d; font-size: 11px;">👉 Click to view WT/ALT/Mismatch details</td></tr>`;

        ["wt", "alt", "wt_mm", "alt_mm"].forEach((key, alleleIdx) => {
            const data = set.alleles?.[key];
            if (!data) return;
            const border = (alleleIdx === 3) ? `border-bottom: 2px solid #adb5bd;` : `border-bottom: 1px solid #eee;`;
            const labelColor = key.includes("mm") ? "#d9534f" : "#0275d8";
            
            html += `<tr class="allele-detail-row ${setRowId}" style="display: none; background-color: #ffffff;">
                <td style="${border}"></td><td style="font-weight: bold; color: ${labelColor}; ${border}">${key.toUpperCase()}</td><td style="border-right: 2px solid #adb5bd; ${border}"></td>
                <td class="mono-cell" style="${border}">${data.oligos.forward.sequence}</td><td class="text-center" style="${border}">${formatNum1(data.oligos.forward.tm)}</td>
                <td class="text-center" style="${border}">${formatNum1(data.oligos.forward.gc)}</td><td class="text-center" style="color:#d9534f;${border}">${formatNum1(data.oligos.forward.hairpin_dg)}</td>
                <td class="text-center" style="border-right: 2px solid #adb5bd;${border} color:#d9534f;">${formatNum1(data.oligos.forward.homodimer_dg)}</td>
                <td class="mono-cell" style="${border}">${data.oligos.reverse.sequence}</td><td class="text-center" style="${border}">${formatNum1(data.oligos.reverse.tm)}</td>
                <td class="text-center" style="${border}">${formatNum1(data.oligos.reverse.gc)}</td><td class="text-center" style="color:#d9534f;${border}">${formatNum1(data.oligos.reverse.hairpin_dg)}</td>
                <td class="text-center" style="border-right: 2px solid #adb5bd;${border} color:#d9534f;">${formatNum1(data.oligos.reverse.homodimer_dg)}</td>
                <td class="text-center" style="border-right: 2px solid #adb5bd;${border} color:#d9534f;">${formatNum1(data.oligos.heterodimer.fr_dg)}</td>
                <td class="text-center" style="${border}">${formatNum1(data.metrics.pair_penalty)}</td><td class="text-center text-nowrap" style="${border}">${data.amplicon_info.genomic_pos}</td>
                <td style="color: ${data.qc_info.is_pass ? '#198754' : '#d9534f'}; font-size: 10px; max-width: 200px; word-wrap: break-word; ${border}">${data.qc_info.fail_reason}</td></tr>`;
        });
    });
    tbody.innerHTML = html;
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
    byId("btn-run-aspcr")?.addEventListener("click", (e) => { e.preventDefault(); runDesign(); });
    
    // 자동 저장
    window.addEventListener("beforeunload", () => {
        TRACKED_INPUTS.forEach(id => {
            const el = byId(id);
            if (el) localStorage.setItem(STORAGE_PREFIX + id, el.type === "checkbox" ? String(el.checked) : el.value);
        });
    });
});