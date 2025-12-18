// app/static/js/toggles.js

document.addEventListener("DOMContentLoaded", () => {
  initAssayTabs();        // qpcr / methyl / as-pcr
  initReferenceToggle();  // hg19 / hg38 ...
  initProbeToggle();      // yes / no + options enable/disable
  initPanelToggles();     // Primer options, QC thresholds 접기/펼치기
  initQcViewToggle();     // QC-only / Design 탭 (있는 페이지에서만)
});

/* =========================
 * 0) Utils
 * ======================= */
function setHidden(el, hidden) {
  if (!el) return;
  el.classList.toggle("hidden", hidden);
}

function setDisabledInside(container, disabled) {
  if (!container) return;
  container.querySelectorAll("input, select, textarea, button").forEach(node => {
    // 버튼까지 disable하면 panel toggle 버튼도 죽을 수 있어서, panel 안쪽만 대상으로 하세요.
    node.disabled = disabled;
  });
}

/* =========================
 * 1) Assay Tabs (primer-tab)
 * - .primer-tab data-assay="qpcr|methyl|as-pcr"
 * - hidden input: #assay-input
 * - form: #design-form (data-action-qpcr/methyl/aspcr 옵션)
 * - panels: #assay-panel-qpcr / #assay-panel-methyl / #assay-panel-aspcr
 * ======================= */
function initAssayTabs() {
  const tabs = document.querySelectorAll(".primer-tab");
  const assayInput = document.getElementById("assay-input");
  const primerTypeInput = document.getElementById("primer-type-input"); // optional
  const form = document.getElementById("design-form");

  if (!tabs.length || !assayInput || !form) return;

  const panelMap = {
    "qpcr": document.getElementById("assay-panel-qpcr"),
    "methyl": document.getElementById("assay-panel-methyl"),
    "as-pcr": document.getElementById("assay-panel-aspcr"),
  };

  function assayToPrimerType(assay) {
    if (assay === "methyl") return "methyl";
    if (assay === "as-pcr") return "as";
    return "default";
  }

  function actionForAssay(assay) {
    if (assay === "methyl") return form.dataset.actionMethyl;
    if (assay === "as-pcr") return form.dataset.actionAspcr;
    return form.dataset.actionQpcr;
  }

  function applyAssay(assay) {
    tabs.forEach(t => t.classList.toggle("active", t.dataset.assay === assay));
    assayInput.value = assay;
    if (primerTypeInput) primerTypeInput.value = assayToPrimerType(assay);

    const nextAction = actionForAssay(assay);
    if (nextAction) form.action = nextAction;

    Object.entries(panelMap).forEach(([k, panel]) => {
      const isActive = (k === assay);
      setHidden(panel, !isActive);
      // 비활성 panel input은 disable 처리 (POST 안되게)
      setDisabledInside(panel, !isActive);
    });

    console.log("[ASSAY]", assay, "action:", form.action);
  }

  applyAssay(assayInput.value || "qpcr");

  tabs.forEach(tab => {
    tab.addEventListener("click", () => {
      const assay = tab.dataset.assay;
      if (!assay) return;
      applyAssay(assay);
    });
  });
}

/* =========================
 * 2) Reference Toggle
 * - hidden input: #reference-input
 * - buttons: .ref-btn-ref data-ref="hg19|hg38"
 * ======================= */
function initReferenceToggle() {
  const referenceInput = document.getElementById("reference-input");
  const refBtns = document.querySelectorAll(".ref-btn-ref");
  if (!referenceInput || !refBtns.length) return;

  function applyRef(ref) {
    referenceInput.value = ref;
    refBtns.forEach(b => b.classList.toggle("active", b.dataset.ref === ref));
    console.log("[REF]", referenceInput.value);
  }

  // ✅ 초기값 기준으로 active 동기화 (서버 렌더 값 반영)
  applyRef(referenceInput.value || "hg38");

  // 클릭 바인딩
  refBtns.forEach(btn => {
    btn.addEventListener("click", () => {
      const ref = btn.dataset.ref;
      if (!ref) return;
      applyRef(ref);
    });
  });
}

/* =========================
 * 3) Probe Toggle (CARD inside)
 * - hidden input: #probe-input (yes/no)
 * - buttons: .probe-btn data-probe="yes|no"
 * - panel: #probe-options (접힘/펼침)
 * - Yes => panel show + panel inputs enabled
 * - No  => panel hide + panel inputs disabled
 * ======================= */
function initProbeToggle() {
  const probeInput = document.getElementById("probe-input");
  const probeBtns = document.querySelectorAll(".probe-btn");
  const probePanel = document.getElementById("probe-options");

  if (!probeBtns.length) return;

  function applyProbeMode(val) {
    if (probeInput) probeInput.value = val;

    probeBtns.forEach(btn => {
      btn.classList.toggle("active", btn.dataset.probe === val);
    });

    // panel show/hide + 내부 입력 enable/disable
    const show = (val === "yes");
    setHidden(probePanel, !show);

    // panel 안쪽 input만 disable (panel null이면 그냥 skip)
    if (probePanel) {
      probePanel.querySelectorAll("input, select, textarea").forEach(node => {
        node.disabled = !show;
      });
    }
  }

  const initial = (probeInput && probeInput.value) ? probeInput.value : "no";
  applyProbeMode(initial);

  probeBtns.forEach(btn => {
    btn.addEventListener("click", () => {
      const val = btn.dataset.probe;
      if (!val) return;
      applyProbeMode(val);
    });
  });
}

/* =========================
 * 4) Panel Toggles (접기/펼치기)
 * - button: [data-toggle="panelId"] data-label="Primer Options"
 * - panel: #panelId
 * - hidden class만 사용 (style.display 섞지 않음)
 * ======================= */
function initPanelToggles() {
  const toggleButtons = document.querySelectorAll("[data-toggle]");

  toggleButtons.forEach(btn => {
    const targetId = btn.dataset.toggle;
    const panel = document.getElementById(targetId);
    if (!panel) return;

    // 초기 텍스트 정합
    const hidden = panel.classList.contains("hidden");
    btn.textContent = `${btn.dataset.label} ${hidden ? "▼" : "▲"}`;

    btn.addEventListener("click", () => {
      panel.classList.toggle("hidden");
      const nowHidden = panel.classList.contains("hidden");
      btn.textContent = `${btn.dataset.label} ${nowHidden ? "▼" : "▲"}`;
    });
  });
}

/* =========================
 * 5) QC View Toggle (Design/QC only tabs)
 * ======================= */
function initQcViewToggle() {
  const tabs = document.querySelectorAll(".qc-tab-btn");
  const designWrapper = document.getElementById("design-wrapper");
  const qcWrapper = document.getElementById("qc-wrapper");

  if (!tabs.length || !designWrapper || !qcWrapper) return;

  function applyQcMode(mode) {
    const isDesign = (mode === "design");
    setHidden(designWrapper, !isDesign);
    setHidden(qcWrapper, isDesign);
    tabs.forEach(t => t.classList.toggle("active", t.dataset.qcMode === mode));
  }

  // 서버에서 기본 active 지정했으면 그걸 따름
  const initial = document.querySelector(".qc-tab-btn.active")?.dataset.qcMode || "design";
  applyQcMode(initial);

  tabs.forEach(tab => {
    tab.addEventListener("click", () => {
      const mode = tab.dataset.qcMode;
      if (!mode) return;
      applyQcMode(mode);
    });
  });
}