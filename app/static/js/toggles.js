// static/js/toggles.js

document.addEventListener("DOMContentLoaded", () => {
  initModeToggle();
  initPrimerTypeToggle();
  initReferenceToggle();
  initProbeToggle();
  initPanelToggles();
  initQcViewToggle();
});

/* =========================
 * 1) Single / Multi Mode 토글
 * ======================= */
function initModeToggle() {
  const modeInput = document.getElementById("mode-input");
  const singleSection = document.getElementById("single-section");
  const multiSection = document.getElementById("multi-section");
  const modeBtns = document.querySelectorAll(".mode-btn");

  // 결과 영역도 있으면 같이 토글
  const singleResults = document.getElementById("single-results");
  const multiResults = document.getElementById("multi-results");

  // 필수 요소 없으면 그냥 스킵
  if (!modeInput || !singleSection || !multiSection || !modeBtns.length) {
    return;
  }

  function applyMode(mode) {
    modeInput.value = mode;

    if (mode === "single") {
      // 폼 영역
      singleSection.style.display = "";
      multiSection.style.display = "none";

      // 결과 영역 (있을 때만)
      if (singleResults) singleResults.style.display = "";
      if (multiResults) multiResults.style.display = "none";
    } else {
      singleSection.style.display = "none";
      multiSection.style.display = "";

      if (singleResults) singleResults.style.display = "none";
      if (multiResults) multiResults.style.display = "";
    }

    // 버튼 active 토글
    modeBtns.forEach(btn => {
      const btnMode = btn.dataset.mode;
      btn.classList.toggle("active", btnMode === mode);
    });

    console.log("[MODE]", modeInput.value);
  }

  // 초기 모드 적용 (서버에서 내려준 값 기반)
  const initialMode = modeInput.value || "single";
  applyMode(initialMode);

  // 클릭 이벤트 바인딩
  modeBtns.forEach(btn => {
    btn.addEventListener("click", () => {
      const mode = btn.dataset.mode;
      if (!mode) return;
      applyMode(mode);
    });
  });
}

/* =========================
 * 2) Primer Type 토글
 * ======================= */
function initPrimerTypeToggle() {
  const primerInput = document.getElementById("primer-type-input");
  const primerBtns = document.querySelectorAll(".primer-btn");

  if (!primerInput || !primerBtns.length) return;

  primerBtns.forEach(btn => {
    btn.addEventListener("click", () => {
      const primerType = btn.getAttribute("data-primer");
      primerInput.value = primerType;

      primerBtns.forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
    });
  });
}

/* =========================
 * 3) Reference 토글
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
 * 4) Probe 토글
 * ======================= */
function initProbeToggle() {
  const probeInput = document.getElementById("probe-input");
  const probeBtns = document.querySelectorAll(".probe-btn");
  const probeOptions = document.getElementById("probe-options");

  if (!probeBtns.length) return;

  // 서버에서 내려준 값이 있으면 사용, 없으면 기본 'no'
  const initialProbe = (probeInput && probeInput.value) ? probeInput.value : "no";

  function applyProbeMode(value) {
    // hidden input 있으면 값 업데이트
    if (probeInput) {
      probeInput.value = value;
    }

    // 버튼 active 상태
    probeBtns.forEach(btn => {
      btn.classList.toggle("active", btn.dataset.probe === value);
    });

    // 옵션 패널 show/hide
    if (probeOptions) {
      if (value === "yes") {
        probeOptions.style.display = "";
      } else {
        probeOptions.style.display = "none";
      }
    }

    console.log("[PROBE]", value);
  }

  // 초기 상태 반영
  applyProbeMode(initialProbe);

  // 클릭 이벤트 바인딩
  probeBtns.forEach(btn => {
    btn.addEventListener("click", () => {
      const val = btn.dataset.probe;
      if (!val) return;
      applyProbeMode(val);
    });
  });
}

/* =========================
 * 5) 옵션 패널 토글 (예: QC Threshold)
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
 * 6) 상단 QC 모드 탭 전환 (QC Only vs Design)
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

  // 초기 모드 (기본 qc_only)
  applyQcMode("design"); // 또는 서버에서 내려준 값 기준으로 바꿀 수 있음

  tabs.forEach(tab => {
    tab.addEventListener("click", () => {
      const mode = tab.dataset.qcMode; // "qc_only" or "design"
      if (!mode) return;
      applyQcMode(mode);
    });
  });
}