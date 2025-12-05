
// static/js/toggles.js
// For browser

document.addEventListener("DOMContentLoaded", () => {
    // === Mode 토글 ===
    const modeInput = document.getElementById("mode-input");
    const singleSection = document.getElementById("single-section");
    const multiSection = document.getElementById("multi-section");
    const modeBtns = document.querySelectorAll(".mode-btn");

    if (modeInput && singleSection && multiSection && modeBtns.length > 0) {
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

    // === Primer Type 토글 ===
    const primerInput = document.getElementById("primer-type-input");
    const primerBtns = document.querySelectorAll(".primer-btn");

    if (primerInput && primerBtns.length > 0) {
        primerBtns.forEach(btn => {
            btn.addEventListener("click", () => {
                const primerType = btn.getAttribute("data-primer");
                primerInput.value = primerType;

                primerBtns.forEach(b => b.classList.remove("active"));
                btn.classList.add("active");
            });
        });
    }

    // === Reference 토글 ===
    const referenceInput = document.getElementById("reference-input");
    const refBtns = document.querySelectorAll(".ref-btn-ref");

    if (referenceInput && refBtns.length > 0) {
        refBtns.forEach(btn => {
            btn.addEventListener("click", () => {
                const ref = btn.getAttribute("data-ref");
                referenceInput.value = ref;

                refBtns.forEach(b => b.classList.remove("active"));
                btn.classList.add("active");
            });
        });
    }

    // === Probe 토글 ===
    const probeInput = document.getElementById("probe-input");
    const probeBtns = document.querySelectorAll(".probe-btn");

    if (probeInput && probeBtns.length > 0) {
        probeBtns.forEach(btn => {
            btn.addEventListener("click", () => {
                const probeVal = btn.getAttribute("data-probe");
                probeInput.value = probeVal;

                probeBtns.forEach(b => b.classList.remove("active"));
                btn.classList.add("active");
            });
        });
    }

    // === QC 패널 토글 ===
    const qcToggleBtn = document.getElementById("qc-toggle-btn");
    const qcPanel = document.getElementById("qc-panel");

    if (qcToggleBtn && qcPanel) {
        qcToggleBtn.addEventListener("click", () => {
            const hidden = qcPanel.classList.toggle("hidden");
            qcToggleBtn.textContent = hidden
                ? "Show QC Thresholds ▼"
                : "Hide QC Thresholds ▲";
        });
    }
});
