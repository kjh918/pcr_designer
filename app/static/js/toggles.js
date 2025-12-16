// app/static/js/toggles.js

document.addEventListener("DOMContentLoaded", () => {
  // initModeToggle();  // ✅ single/multi 제거했으면 호출 안 해도 됨
  initPrimerTypeToggle();
  initReferenceToggle();
  initProbeToggle();
  initPanelToggles();
  initQcViewToggle();
});

function initPrimerTypeToggle() {
  const tabs = document.querySelectorAll(".primer-tab");
  const assayInput = document.getElementById("assay-input");
  const primerTypeInput = document.getElementById("primer-type-input"); // optional
  const form = document.getElementById("design-form");

  if (!tabs.length || !assayInput || !form) return;

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

  function setVisibilityAndDisable(activeAssay) {
    const groups = [
      { assay: "qpcr",  ids: ["assay-input-qpcr",  "assay-panel-qpcr"] },
      { assay: "methyl", ids: ["assay-input-methyl", "assay-panel-methyl"] },
      { assay: "as-pcr", ids: ["assay-input-aspcr", "assay-panel-aspcr"] },
    ];

    groups.forEach(g => {
      const isActive = g.assay === activeAssay;

      g.ids.forEach(id => {
        const el = document.getElementById(id);
        if (!el) return;

        el.style.display = isActive ? "" : "none";
        el.classList.toggle("hidden", !isActive);

        el.querySelectorAll("input, select, textarea").forEach(node => {
          node.disabled = !isActive;
        });
      });
    });
  }

  function applyAssay(assay) {
    // 1) active 탭 표시
    tabs.forEach(t => t.classList.toggle("active", t.dataset.assay === assay));

    // 2) hidden 업데이트
    assayInput.value = assay;
    if (primerTypeInput) primerTypeInput.value = assayToPrimerType(assay);

    // 3) form action 변경
    const nextAction = actionForAssay(assay);
    if (nextAction) form.action = nextAction;

    // 4) 패널 전환
    setVisibilityAndDisable(assay);

    console.log("[ASSAY]", assay, "action:", form.action);
  }

  // 초기 적용(서버 렌더 hidden 값 기준)
  applyAssay(assayInput.value || "qpcr");

  // 클릭 바인딩
  tabs.forEach(tab => {
    tab.addEventListener("click", () => {
      const assay = tab.dataset.assay;
      if (!assay) return;
      applyAssay(assay);
    });
  });
}

  function applyAssayFromPrimer(primerType) {
    const assay = primerToAssay(primerType);

    // hidden 업데이트
    assayInput.value = assay;
    if (primerInput) primerInput.value = primerType;

    // 버튼 active
    primerBtns.forEach(b => b.classList.toggle("active", b.dataset.primer === primerType));

    // form action 변경
    const nextAction = actionForAssay(assay);
    if (nextAction) form.action = nextAction;

    // 패널 전환
    setVisibilityAndDisable(assay);

    console.log("[TAB]", primerType, "=>", assay, "action:", form.action);
  }

  // 초기 상태 결정: hidden assay 기준으로 primerType 역매핑
  function assayToPrimer(assay) {
    if (assay === "methyl") return "methyl";
    if (assay === "as-pcr") return "as";
    return "default";
  }

  const initialAssay = assayInput.value || "qpcr";
  const initialPrimerType =
    (primerInput && primerInput.value) ||
    (document.querySelector(".primer-btn.active")?.dataset.primer) ||
    assayToPrimer(initialAssay);

  applyAssayFromPrimer(initialPrimerType);

  // 클릭 이벤트
  primerBtns.forEach(btn => {
    btn.addEventListener("click", () => {
      const primerType = btn.dataset.primer;
      if (!primerType) return;
      applyAssayFromPrimer(primerType);
    });
  });
}

/* =========================
 * 2) Reference 토글
 * ======================= */
function initReferenceToggle() {
  const referenceInput = document.getElementById("reference-input");
  const refBtns = document.querySelectorAll(".ref-btn-ref");

  if (!referenceInput || !refBtns.length) return;

  refBtns.forEach(btn => {
    btn.addEventListener("click", () => {
      const ref = btn.getAttribute("data-ref");
      referenceInput.value = ref;

      refBtns.forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
    });
  });
}

/* =========================
 * 3) Probe 토글
 * ======================= */
function initProbeToggle() {
  const probeInput = document.getElementById("probe-input");
  const probeBtns = document.querySelectorAll(".probe-btn");
  const probeOptions = document.getElementById("probe-options");

  if (!probeBtns.length) return;

  const initialProbe = (probeInput && probeInput.value) ? probeInput.value : "no";

  function applyProbeMode(value) {
    if (probeInput) probeInput.value = value;

    probeBtns.forEach(btn => {
      btn.classList.toggle("active", btn.dataset.probe === value);
    });

    if (probeOptions) {
      probeOptions.style.display = (value === "yes") ? "" : "none";
    }
  }

  applyProbeMode(initialProbe);

  probeBtns.forEach(btn => {
    btn.addEventListener("click", () => {
      const val = btn.dataset.probe;
      if (!val) return;
      applyProbeMode(val);
    });
  });
}

/* =========================
 * 4) 옵션 패널 토글
 * ======================= */
function initPanelToggles() {
  const toggleButtons = document.querySelectorAll("[data-toggle]");

  toggleButtons.forEach(btn => {
    const targetId = btn.getAttribute("data-toggle");
    const panel = document.getElementById(targetId);
    if (!panel) return;

    btn.addEventListener("click", () => {
      const isHidden =
        panel.style.display === "none" ||
        (panel.style.display === "" && panel.classList.contains("hidden"));

      if (isHidden) {
        panel.style.display = "";
        btn.textContent = `${btn.dataset.label} ▲`;
      } else {
        panel.style.display = "none";
        btn.textContent = `${btn.dataset.label} ▼`;
      }
    });
  });
}

/* =========================
 * 5) QC 모드 탭 전환 (있을 때만)
 * ======================= */
function initQcViewToggle() {
  const tabs = document.querySelectorAll(".qc-tab-btn");
  const designWrapper = document.getElementById("design-wrapper");
  const qcWrapper = document.getElementById("qc-wrapper");

  if (!tabs.length || !designWrapper || !qcWrapper) return;

  function applyQcMode(mode) {
    if (mode === "design") {
      designWrapper.style.display = "";
      qcWrapper.style.display = "none";
    } else {
      designWrapper.style.display = "none";
      qcWrapper.style.display = "";
    }

    tabs.forEach(t => {
      t.classList.toggle("active", t.dataset.qcMode === mode);
    });
  }

  applyQcMode("design");

  tabs.forEach(tab => {
    tab.addEventListener("click", () => {
      const mode = tab.dataset.qcMode;
      if (!mode) return;
      applyQcMode(mode);
    });
  });
}
