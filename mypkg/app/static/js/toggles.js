// static/js/toggles.js

document.addEventListener("DOMContentLoaded", () => {
    initModeToggle();
    initPrimerTypeToggle();
    initReferenceToggle();
    initProbeToggle();
    initQcPanelToggle();
    initQcViewToggle();   // ← 이 줄이 실제로 있는지!
});

// === 1) Single / Multi 설계 모드 토글 ===
function initModeToggle() {
    const modeInput = document.getElementById("mode-input");
    const singleSection = document.getElementById("single-section");
    const multiSection = document.getElementById("multi-section");
    const modeBtns = document.querySelectorAll(".mode-btn");

    if (!modeInput || !singleSection || !multiSection || !modeBtns.length) return;

    modeBtns.forEach(btn => {
        btn.addEventListener("click", () => {
            const mode = btn.getAttribute("data-mode");
            modeInput.value = mode;

            modeBtns.forEach(b => b.classList.remove("active"));
            btn.classList.add("active");

            if (mode === "single") {
                singleSection.classList.remove("hidden");
                multiSection.classList.add("hidden");
            } else {
                singleSection.classList.add("hidden");
                multiSection.classList.remove("hidden");
            }
        });
    });
}

// === 2) Primer Type 토글 ===
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

// === 3) Reference 토글 ===
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

// === 4) Probe 토글 ===
function initProbeToggle() {
    const probeInput = document.getElementById("probe-input");
    const probeBtns = document.querySelectorAll(".probe-btn");

    if (!probeInput || !probeBtns.length) return;

    probeBtns.forEach(btn => {
        btn.addEventListener("click", () => {
            const probeVal = btn.getAttribute("data-probe");
            probeInput.value = probeVal;

            probeBtns.forEach(b => b.classList.remove("active"));
            btn.classList.add("active");
        });
    });
}

// === 5) QC Threshold 패널 토글 ===
function initQcPanelToggle() {
    const qcToggleBtn = document.getElementById("qc-toggle-btn");
    const qcPanel = document.getElementById("qc-panel");

    if (!qcToggleBtn || !qcPanel) return;

    qcToggleBtn.addEventListener("click", () => {
        const hidden = qcPanel.classList.toggle("hidden");
        qcToggleBtn.textContent = hidden
            ? "Show QC Thresholds ▼"
            : "Hide QC Thresholds ▲";
    });
}
// === 6) Primer Design / QC Only wrapper 전환 ===
function initQcViewToggle() {
    const tabs = document.querySelectorAll(".qc-tab-btn");
    const designWrapper = document.getElementById("design-wrapper");
    const qcWrapper = document.getElementById("qc-wrapper");

    if (!tabs.length || !designWrapper || !qcWrapper) return;

    tabs.forEach(tab => {
        tab.addEventListener("click", () => {
            const mode = tab.dataset.qcMode;  // "qc_only" or "design"

            tabs.forEach(t => t.classList.remove("active"));
            tab.classList.add("active");

            if (mode === "design") {
                designWrapper.classList.remove("hidden");
                qcWrapper.classList.add("hidden");
            } else {
                designWrapper.classList.add("hidden");
                qcWrapper.classList.remove("hidden");
            }
        });
    });
}
