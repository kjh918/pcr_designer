// === Mode 토글 ===
const modeInput = document.getElementById("mode-input");
const singleSection = document.getElementById("single-section");
const multiSection = document.getElementById("multi-section");
const modeBtns = document.querySelectorAll(".mode-btn");

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

// === Primer Type 토글 ===
const primerInput = document.getElementById("primer-type-input");
const primerBtns = document.querySelectorAll(".primer-btn");

primerBtns.forEach(btn => {
    btn.addEventListener("click", () => {
        const primerType = btn.getAttribute("data-primer");
        primerInput.value = primerType;

        primerBtns.forEach(b => b.classList.remove("active"));
        btn.classList.add("active");
    });
});

// === Reference 토글 ===
const referenceInput = document.getElementById("reference-input");
const refBtns = document.querySelectorAll(".ref-btn-ref");

refBtns.forEach(btn => {
    btn.addEventListener("click", () => {
        const ref = btn.getAttribute("data-ref");
        referenceInput.value = ref;

        refBtns.forEach(b => b.classList.remove("active"));
        btn.classList.add("active");
    });
});


// === Probe 토글 ===
const probeInput = document.getElementById("probe-input");
const probeBtns = document.querySelectorAll(".probe-btn");

probeBtns.forEach(btn => {
    btn.addEventListener("click", () => {
        const probeVal = btn.getAttribute("data-probe");
        probeInput.value = probeVal;

        probeBtns.forEach(b => b.classList.remove("active"));
        btn.classList.add("active");
    });
});
