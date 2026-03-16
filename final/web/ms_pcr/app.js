console.log("🚀 [v2.0] MS-PCR app.js 파일이 성공적으로 로드되었습니다! (캐시 갱신 완료)");

const state = { lastResults: null, lastPayload: null };
const STORAGE_PREFIX = "mspcr_";
// 🔥 프론트엔드 요청 주소 (백엔드와 정확히 일치해야 합니다)
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
  if (!el) {
      console.warn(`⚠️ [경고] HTML에서 ID가 '${id}'인 요소를 찾을 수 없습니다. 기본값(${def})을 사용합니다.`);
      return def;
  }
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
  console.log("🛠️ 페이로드 조립 시작...");
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
  console.log("🎨 결과 렌더링 시작...");
  const tbody = byId("candidate-tbody");
  if (!tbody) {
      console.error("❌ 'candidate-tbody' 테이블을 찾을 수 없습니다.");
      return;
  }
  tbody.innerHTML = "";

  const listData = results?.single_filtered_amplicons || [];
  if (!listData.length) {
    tbody.innerHTML = '<tr><td colspan="10" style="text-align:center; padding:40px;">조건을 만족하는 MS-PCR 세트를 찾지 못했습니다.</td></tr>';
    return;
  }
listData.forEach((set, idx) => {
    const rank = set.rank ?? (idx + 1);
    const m = set.alleles?.M;
    const u = set.alleles?.U;
    if (!m || !u) return;

    // 🔥 QC 통과 여부 뱃지 생성 (마우스 올리면 실패 사유 표시)
    const mQcBadge = m.is_qc_pass ? '<span style="color: #28a745; font-weight: 900;">PASS</span>' : `<span style="color: #dc3545; font-weight: 900;" title="${m.qc_log || ''}">FAIL</span>`;
    const uQcBadge = u.is_qc_pass ? '<span style="color: #28a745; font-weight: 900;">PASS</span>' : `<span style="color: #dc3545; font-weight: 900;" title="${u.qc_log || ''}">FAIL</span>`;

    // M-Allele
    const trM = document.createElement("tr");
    trM.style.cursor = "pointer";
    trM.style.backgroundColor = "#fff5f5"; // 연한 빨강 배경만 남김
    trM.innerHTML = `
      <td class="text-center" rowspan="2" style="vertical-align: middle; border-bottom: 2px solid #ccc; border-right: 1px solid #ddd;"><b>${rank}</b></td>
      <td class="text-center"><span style="background:#d9534f; color:white; padding:2px 6px; border-radius:3px; font-size:11px; font-weight:bold;">M-Allele</span></td>
      <td class="mono-cell">${m.oligos.forward.sequence}</td>
      <td class="mono-cell">${m.oligos.reverse.sequence}</td>
      <td class="text-center">${formatNum1(m.oligos.forward.tm)}</td>
      <td class="text-center">${formatNum1(m.oligos.reverse.tm)}</td>
      <td class="text-center">${formatNum1(m.oligos.forward.gc)}%</td>
      <td class="text-center">${formatNum1(m.oligos.reverse.gc)}%</td>
      <td class="text-center"><b>${m.amplicon_info.sequence.length} bp</b></td>
      <td class="text-center">${formatNum1(m.metrics.pair_penalty)}</td>
      <td class="text-center">${mQcBadge}</td>
    `;

    // U-Allele
    const trU = document.createElement("tr");
    trU.style.cursor = "pointer";
    trU.style.backgroundColor = "#f0f7ff"; // 연한 파랑 배경만 남김
    trU.style.borderBottom = "2px solid #ccc";
    trU.innerHTML = `
      <td class="text-center"><span style="background:#0275d8; color:white; padding:2px 6px; border-radius:3px; font-size:11px; font-weight:bold;">U-Allele</span></td>
      <td class="mono-cell">${u.oligos.forward.sequence}</td>
      <td class="mono-cell">${u.oligos.reverse.sequence}</td>
      <td class="text-center">${formatNum1(u.oligos.forward.tm)}</td>
      <td class="text-center">${formatNum1(u.oligos.reverse.tm)}</td>
      <td class="text-center">${formatNum1(u.oligos.forward.gc)}%</td>
      <td class="text-center">${formatNum1(u.oligos.reverse.gc)}%</td>
      <td class="text-center"><b>${u.amplicon_info.sequence.length} bp</b></td>
      <td class="text-center">${formatNum1(u.metrics.pair_penalty)}</td>
      <td class="text-center">${uQcBadge}</td>
    `;

    const onClickRow = () => {
      document.querySelectorAll("#candidate-tbody tr").forEach((r) => r.classList.remove("selected-row"));
      trM.classList.add("selected-row");
      trU.classList.add("selected-row");
      const alnView = byId("alignment-view");
      if(alnView) alnView.innerText = m.alignment_text_block || m.alignment_visual?.join("\n") || "No alignment data.";
    };

    trM.addEventListener("click", onClickRow);
    trU.addEventListener("click", onClickRow);

    tbody.appendChild(trM);
    tbody.appendChild(trU);
  });

  
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
  console.log("▶️ DESIGN START 버튼 클릭됨!");

  const runBtn = byId("btn-run-mspcr");
  if (!runBtn) {
      alert("HTML에서 'btn-run-mspcr' 버튼을 찾을 수 없습니다! HTML 코드를 확인해주세요.");
      return;
  }

  saveInputsToStorage();
  const payload = buildPayload();
  console.log("📦 전송할 데이터:", payload);
  
  const err = validatePayload(payload);
  if (err) {
      console.error("❌ 검증 실패:", err);
      alert(err);
      return;
  }

  const originalText = runBtn.innerText;
  runBtn.innerText = "⏳ RUNNING...";
  runBtn.disabled = true;

  try {
    console.log(`🌐 API 호출 중... (${API_URL})`);
    const res = await fetch(API_URL, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });

    console.log("🌐 서버 응답 코드:", res.status);
    if (!res.ok) {
        const errorJson = await res.json().catch(()=>({}));
        throw new Error(`서버 에러 (${res.status}): ${errorJson.detail || "원인 불명"}`);
    }
    
    const results = await res.json();
    console.log("📥 백엔드 응답 데이터:", results);
    
    if (results.status === "fail") throw new Error(results.error || results.reason || "Unknown Error");

    state.lastPayload = payload;
    state.lastResults = results;

    renderResults(results);
    renderMetadata(results, payload);
    renderSummary(payload, results);
    
    console.log("✅ 디자인 완료 처리 성공!");

  } catch (e) {
    console.error("🔥 Fetch 에러 발생:", e);
    alert("요청 실패: " + e.message);
    const statusEl = byId("summary-status");
    if(statusEl) {
        statusEl.innerText = "ERROR";
        statusEl.style.color = "red";
    }
  } finally {
    runBtn.innerText = originalText;
    runBtn.disabled = false;
    console.log("⏹️ 실행 흐름 종료. 버튼 활성화.");
  }
}

document.addEventListener("DOMContentLoaded", () => {
  console.log("✅ DOMContentLoaded 이벤트 발생! 요소들을 연결합니다.");
  restoreInputsFromStorage();

  const btnRun = byId("btn-run-mspcr");
  if (btnRun) {
      console.log("✅ 실행 버튼 발견! 이벤트를 연결합니다.");
      btnRun.addEventListener("click", (e) => { 
          e.preventDefault(); 
          runDesign(); 
      });
  } else {
      console.error("❌ 'btn-run-mspcr' 버튼을 화면에서 찾을 수 없습니다! HTML 코드를 다시 확인하세요.");
      alert("에러: 디자인 시작 버튼을 화면에서 찾을 수 없습니다.");
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
      console.log("✅ 튜토리얼 데이터가 로드되었습니다.");
  });

  console.log("🎉 모든 초기화 완료.");
});