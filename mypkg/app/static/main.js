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

const btn = document.getElementById("qc-toggle-btn");
const panel = document.getElementById("qc-panel");

btn.addEventListener("click", () => {
    const hidden = panel.classList.toggle("hidden");
    btn.textContent = hidden ? "Show QC Thresholds ▼" : "Hide QC Thresholds ▲";
})

const express = require('express');
const path = require('path');

const app = express();
const PORT = process.env.PORT || 3000;

// body parser
app.use(express.json());
app.use(express.urlencoded({ extended: true }));

// 정적 파일 (style.css, 기존 html 등)
app.use(express.static(path.join(__dirname, 'public')));

// QC 설정 받는 엔드포인트 (BLAST 관련 없음)
