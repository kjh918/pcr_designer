console.log("🚀 [v4.0] MS-PCR app.js 로드 완료! (스크롤 및 QC 사유 직접 표기 반영)");

const state = { lastResults: null, lastPayload: null };
const STORAGE_PREFIX = "mspcr_";
const API_URL = "http://192.168.0.35:9000/api/design/mspcr"; 

const TRACKED_INPUTS = [
  "mspcr-design-name", "mspcr-raw-sequence", "input-ref-genome",
  "min_amplicon_length", "max_amplicon_length",
  "primer_min_length", "primer_opt_length", "primer_max_length",
  "primer_min_tm", "primer_opt_tm", "primer_max_tm",
  "primer_min_gc", "primer_opt_gc", "primer_max_gc",
  "qc_hairpin_min_dg", "qc_homodimer_min_dg", "qc_heterodimer_min_dg",
  "qc_min_identity", "qc_min_hit_length", "qc_blast_max_alignments",
  "qc_use_ispcr_check", "qc_primer_max_diff_tm",
];

const $ = (sel) => document.querySelector(sel);
function byId(id) { return document.getElementById(id); }

function formatNum1(val) {
  if (val === "-" || val == null || Number.isNaN(val)) return "-";
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
    design_name: getVal("mspcr-design-name", "string", "MS-PCR_Project"),
    sequence: getVal("mspcr-raw-sequence", "string", "").replace(/\s/g, "").toUpperCase(),
    reference: getVal("input-ref-genome", "string", "hg38"),
    top_k: 5,
    min_amplicon_length: getVal("min_amplicon_length", "int", 80),
    max_amplicon_length: getVal("max_amplicon_length", "int", 200),
    primer_min_length: getVal("primer_min_length", "int", 18),
    primer_opt_length: getVal("primer_opt_length", "int", 22),
    primer_max_length: getVal("primer_max_length", "int", 30),
    primer_min_tm: getVal("primer_min_tm", "float", 52.0),
    primer_opt_tm: getVal("primer_opt_tm", "float", 58.0),
    primer_max_tm: getVal("primer_max_tm", "float", 65.0),
    primer_min_gc: getVal("primer_min_gc", "float", 30.0),
    primer_opt_gc: getVal("primer_opt_gc", "float", 50.0),
    primer_max_gc: getVal("primer_max_gc", "float", 70.0),
    qc_hairpin_min_dg: getVal("qc_hairpin_min_dg", "float", -5.0),
    qc_homodimer_min_dg: getVal("qc_homodimer_min_dg", "float", -6.0),
    qc_heterodimer_min_dg: getVal("qc_heterodimer_min_dg", "float", -6.0),
    qc_min_identity: getVal("qc_min_identity", "float", 90.0),
    qc_min_hit_length: getVal("qc_min_hit_length", "int", 13),
    qc_blast_max_alignments: getVal("qc_blast_max_alignments", "int", 50),
    qc_use_ispcr_check: getVal("qc_use_ispcr_check", "bool", false),
    qc_primer_max_diff_tm: getVal("qc_primer_max_diff_tm", "float", 3.0)
  };
}

function validatePayload(p) {
  if (!p.sequence || p.sequence.length < 20) return "올바른 서열을 입력해주세요 (최소 20bp 이상).";
  if (!p.sequence.includes("[") || !p.sequence.includes("]")) return "타겟 CpG를 대괄호 [CG]로 감싸주세요!\n(예: ATGCATGC[CG]ATGCAT)";
  return null;
}

function renderResults(results) {
  const tbody = byId("candidate-tbody");
  if (!tbody) return;
  tbody.innerHTML = "";

  let listData = results?.single_filtered_amplicons || [];
  if (!listData.length) {
    tbody.innerHTML = '<tr><td colspan="13" style="text-align:center; padding:40px;">조건을 만족하는 MS-PCR 세트를 찾지 못했습니다.</td></tr>';
    return;
  }

  // 🔥 1. 정렬: PASS인 세트가 무조건 위로 오도록 오름차순(우선순위) 정렬
  listData.sort((a, b) => {
      const aPass = a.set_qc_pass ? 1 : 0;
      const bPass = b.set_qc_pass ? 1 : 0;
      if (aPass !== bPass) return bPass - aPass; // PASS(1)가 FAIL(0)보다 위로
      return (a.rank || 0) - (b.rank || 0);      // 같은 상태면 원래 랭크(Penalty) 순 유지
  });

  // 🔥 2. 데이터 수 제한: 최대 100개 세트 (화면 부하 방지)
  listData = listData.slice(0, 100);

  listData.forEach((set, idx) => {
    const rank = idx + 1; // 재정렬되었으므로 1부터 새로 넘버링
    const m = set.alleles?.M;
    const u = set.alleles?.U;
    if (!m || !u) return;

    // 🔥 3. QC 사유 직접 추출 및 렌더링
    const mPass = m.qc_info?.is_pass;
    const mFailReason = m.qc_info?.fail_reason || "Unknown Fail";
    const uPass = u.qc_info?.is_pass;
    const uFailReason = u.qc_info?.fail_reason || "Unknown Fail";

    const mQcBadge = mPass 
        ? '<span style="color: #28a745; font-weight: 900;">PASS</span>' 
        : `<span style="color: #dc3545; font-weight: 900;">FAIL</span><div style="font-size:10px; color:#888; line-height:1.2; margin-top:2px; word-break:keep-all;">${mFailReason}</div>`;
    
    const uQcBadge = uPass 
        ? '<span style="color: #28a745; font-weight: 900;">PASS</span>' 
        : `<span style="color: #dc3545; font-weight: 900;">FAIL</span><div style="font-size:10px; color:#888; line-height:1.2; margin-top:2px; word-break:keep-all;">${uFailReason}</div>`;

    // M-Allele 생성
    const trM = document.createElement("tr");
    trM.style.cursor = "pointer";
    trM.style.backgroundColor = "#fff5f5";
    trM.innerHTML = `
      <td class="text-center" rowspan="2" style="vertical-align: middle; border-bottom: 2px solid #ccc; border-right: 1px solid #ddd; background: #fff;"><b>${rank}</b></td>
      <td class="text-center" style="border-bottom: 1px solid #f5c6cb;"><span style="background:#d9534f; color:white; padding:2px 6px; border-radius:3px; font-size:11px; font-weight:bold;">M-Allele</span></td>
      <td class="mono-cell" style="border-bottom: 1px solid #f5c6cb;">${m.oligos.forward.sequence}</td>
      <td class="mono-cell" style="border-bottom: 1px solid #f5c6cb;">${m.oligos.reverse.sequence}</td>
      <td class="text-center" style="border-bottom: 1px solid #f5c6cb;">${formatNum1(m.oligos.forward.tm)}</td>
      <td class="text-center" style="border-bottom: 1px solid #f5c6cb;">${formatNum1(m.oligos.reverse.tm)}</td>
      <td class="text-center" style="border-bottom: 1px solid #f5c6cb;">${formatNum1(m.oligos.forward.gc)}%</td>
      <td class="text-center" style="border-bottom: 1px solid #f5c6cb;">${formatNum1(m.oligos.reverse.gc)}%</td>
      <td class="text-center" style="border-bottom: 1px solid #f5c6cb;"><b>${m.amplicon_info.sequence.length} bp</b></td>
      <td class="text-center" style="border-bottom: 1px solid #f5c6cb; color:#1e707a; font-weight:bold;">${formatNum1(m.amplicon_info.tm)}</td>
      <td class="text-center" style="border-bottom: 1px solid #f5c6cb; color:#1e707a;">${formatNum1(m.amplicon_info.gc)}%</td>
      <td class="text-center" style="border-bottom: 1px solid #f5c6cb; color:#1e707a;">${formatNum1(m.amplicon_info.cpg_count)}</td>
      <td class="text-center" style="border-bottom: 1px solid #f5c6cb;">${mQcBadge}</td>
    `;

    // U-Allele 생성
    const trU = document.createElement("tr");
    trU.style.cursor = "pointer";
    trU.style.backgroundColor = "#f0f7ff";
    trU.innerHTML = `
      <td class="text-center" style="border-bottom: 2px solid #ccc;"><span style="background:#0275d8; color:white; padding:2px 6px; border-radius:3px; font-size:11px; font-weight:bold;">U-Allele</span></td>
      <td class="mono-cell" style="border-bottom: 2px solid #ccc;">${u.oligos.forward.sequence}</td>
      <td class="mono-cell" style="border-bottom: 2px solid #ccc;">${u.oligos.reverse.sequence}</td>
      <td class="text-center" style="border-bottom: 2px solid #ccc;">${formatNum1(u.oligos.forward.tm)}</td>
      <td class="text-center" style="border-bottom: 2px solid #ccc;">${formatNum1(u.oligos.reverse.tm)}</td>
      <td class="text-center" style="border-bottom: 2px solid #ccc;">${formatNum1(u.oligos.forward.gc)}%</td>
      <td class="text-center" style="border-bottom: 2px solid #ccc;">${formatNum1(u.oligos.reverse.gc)}%</td>
      <td class="text-center" style="border-bottom: 2px solid #ccc;"><b>${u.amplicon_info.sequence.length} bp</b></td>
      <td class="text-center" style="border-bottom: 2px solid #ccc; color:#1e707a; font-weight:bold;">${formatNum1(u.amplicon_info.tm)}</td>
      <td class="text-center" style="border-bottom: 2px solid #ccc; color:#1e707a;">${formatNum1(u.amplicon_info.gc)}%</td>
      <td class="text-center" style="border-bottom: 2px solid #ccc; color:#1e707a;">${formatNum1(u.amplicon_info.cpg_count)}</td>
      <td class="text-center" style="border-bottom: 2px solid #ccc;">${uQcBadge}</td>
    `;

    const onClickRow = () => {
      document.querySelectorAll("#candidate-tbody tr").forEach((r) => r.classList.remove("selected-row"));
      trM.classList.add("selected-row");
      trU.classList.add("selected-row");
      const alnView = byId("alignment-view");
      if(alnView) alnView.innerText = m.alignment_text_block || "No alignment data.";
    };

    trM.addEventListener("click", onClickRow);
    trU.addEventListener("click", onClickRow);

    tbody.appendChild(trM);
    tbody.appendChild(trU);
  });

  document.querySelector("#candidate-tbody tr")?.click();
}

function renderMetadata(results, payload) {
  const panel = byId("metadata-panel");
  const container = byId("metadata-table");
  if (!panel || !container || !results || !payload) return;

  panel.style.display = "block";
  const filtered = results.export_meta?.filtered_count || 0;
  const conv = results.conversion_info || {};

  container.innerHTML = `
    <div class="meta-report-grid">
      <div class="meta-card">
        <div class="card-header">📍 TARGET INFO</div>
        <div class="card-item"><span>Target Pos</span> <b class="text-danger">${results.target_info?.target_range || "-"}</b></div>
        <div class="card-item"><span>Total CpG</span> <b>${conv.total_cpg_count || 0} sites</b></div>
        <div class="card-item"><span>Amp Size</span> <b>${payload.min_amplicon_length}-${payload.max_amplicon_length} bp</b></div>
        <div class="card-item"><span>Generated Sets</span> <b>${filtered} Sets</b></div>
      </div>
      <div class="meta-card">
        <div class="card-header">🔄 CONVERSION RULE</div>
        <div class="card-item"><span>M-Allele</span> <b>Target CpG 유지, 나머지 C->T</b></div>
        <div class="card-item"><span>U-Allele</span> <b>모든 C를 T로 변환</b></div>
        <div class="card-item"><span>3' Anchoring</span> <b style="color: #28a745;">APPLIED</b></div>
      </div>
      <div class="meta-card">
        <div class="card-header">🛡️ BLAST & SPECIFICITY</div>
        <div class="card-item"><span>Min Identity</span> <b>${payload.qc_min_identity}%</b></div>
        <div class="card-item"><span>Max Align</span> <b>${payload.qc_blast_max_alignments} Hits</b></div>
        <div class="card-item"><span>Max ΔTm</span> <b>${payload.qc_primer_max_diff_tm} ℃</b></div>
      </div>
    </div>
  `;
}

function renderSummary(payload, results) {
  const now = new Date();
  const summaryDate = byId("summary-date");
  if(summaryDate) summaryDate.innerText = now.toLocaleDateString() + " " + now.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });
  
  const summaryMut = byId("summary-mutation");
  if(summaryMut) summaryMut.innerText = `[CG] Target`; 

  const statusEl = byId("summary-status");
  if(!statusEl) return;
  
  if (results?.status === "success") { statusEl.innerText = "COMPLETED"; statusEl.style.color = "green"; }
  else if (results?.status) { statusEl.innerText = "FAILED"; statusEl.style.color = "red"; }
}

async function runDesign() {
  const runBtn = byId("btn-run-mspcr");
  if (!runBtn) return;

  saveInputsToStorage();
  const payload = buildPayload();
  const err = validatePayload(payload);
  if (err) { alert(err); return; }

  const originalText = runBtn.innerText;
  runBtn.innerText = "⏳ RUNNING...";
  runBtn.disabled = true;

  try {
    const res = await fetch(API_URL, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });

    if (!res.ok) {
        const errorJson = await res.json().catch(()=>({}));
        throw new Error(`서버 에러 (${res.status}): ${errorJson.detail || "원인 불명"}`);
    }
    
    const results = await res.json();
    if (results.status === "fail") throw new Error(results.error || results.reason || "Unknown Error");

    state.lastPayload = payload;
    state.lastResults = results;

    renderResults(results);
    renderMetadata(results, payload);
    renderSummary(payload, results);

  } catch (e) {
    alert("요청 실패: " + e.message);
    const statusEl = byId("summary-status");
    if(statusEl) {
        statusEl.innerText = "ERROR";
        statusEl.style.color = "red";
    }
  } finally {
    runBtn.innerText = originalText;
    runBtn.disabled = false;
  }
}

document.addEventListener("DOMContentLoaded", () => {
  restoreInputsFromStorage();

  const btnRun = byId("btn-run-mspcr");
  if (btnRun) {
      btnRun.addEventListener("click", (e) => { 
          e.preventDefault(); 
          runDesign(); 
      });
  }

  byId("btn-load-tutorial")?.addEventListener("click", (e) => {
      e.preventDefault();
      const txtName = byId("mspcr-design-name");
      const txtSeq = byId("mspcr-raw-sequence");
      if(txtName) txtName.value = "MGMT_Promoter_Tutorial";
      if(txtSeq) txtSeq.value = "ACTGCTAGCTGATCGATCGATCGACTGAC[CG]TCGATCGATCGACTGCATCGATCGATCGACTGCTAGCTGATCGATCGATCGACTAGCTAGCTAGCTAGCGCGATCGACTAGCTGATC";
      
      const btn = e.target;
      btn.innerText = "✅ Loaded!";
      setTimeout(() => btn.innerText = "🧪 Tutorial", 2000);
  });
});