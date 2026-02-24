/**
 * qpcr_app.js
 * - 입력값 자동 저장/복구 (LocalStorage)
 * - FastAPI 서버와 통신 (Real Fetch)
 * - 결과 테이블 렌더링
 */

// 1. 추적할 입력창 ID 목록
const TRACKED_INPUTS = [
    // --- Neccessary Inputs ---
    "input-ref-genome","input-chrom", "input-start", "input-end",
    "input-ref-base", "input-alt-base", "input-strand",
    
    // --- Primer Options ---
    "min_amplicon_length", "max_amplicon_length",
    "primer_min_length", "primer_opt_length", "primer_max_length",
    "primer_min_tm", "primer_opt_tm", "primer_max_tm",
    "primer_min_gc", "primer_opt_gc", "primer_max_gc",

    // --- Probe Options ---
    "probe_min_length", "probe_opt_length", "probe_max_length",
    "min_primer_probe_tm_diff", "max_primer_probe_tm_diff",
    "probe_min_tm", "probe_opt_tm", "probe_max_tm",
    "probe_min_gc", "probe_opt_gc", "probe_max_gc"
    // 필요시 여기에 다른 ID 추가 (예: "min_amplicon_length" 등)
];

// 2. 브라우저 저장소에 값 저장
function saveInputsToStorage() {
    TRACKED_INPUTS.forEach(id => {
        const inputEl = document.getElementById(id);
        if (inputEl && inputEl.value) {
            localStorage.setItem("qpcr_" + id, inputEl.value);
        }
    });
    console.log("💾 Input values saved to LocalStorage.");
}

// 3. 브라우저 저장소에서 값 복구
function restoreInputsFromStorage() {
    TRACKED_INPUTS.forEach(id => {
        const savedValue = localStorage.getItem("qpcr_" + id);
        const inputEl = document.getElementById(id);
        if (savedValue !== null && inputEl) {
            inputEl.value = savedValue;
        }
    });
    console.log("📂 Input values restored from LocalStorage.");
}
/**
 * 서버에서 받은 JSON 데이터를 결과 테이블에 매핑합니다.
 * (수정됨: 데이터가 리스트가 아닐 경우를 대비한 안전장치 추가)
 */
function renderResults(data) {
    const tbody = document.getElementById("candidate-tbody");
    if (!tbody) return;
    tbody.innerHTML = ""; 

    console.log("📦 [Web] 렌더링할 데이터 원본:", data); // F12 콘솔에서 확인용

    // 💡 핵심 수정: 데이터가 배열(List)이 아니면 배열로 변환하거나 내부를 찾습니다.
    let listData = [];
    
    if (Array.isArray(data)) {
        // 1. 정상적인 리스트 형태 [...]
        listData = data;
    } else if (data && typeof data === 'object') {
        // 2. 딕셔너리 형태 {...} 라면, 혹시 내부에 리스트가 숨어있나 확인
        if (data.results && Array.isArray(data.results)) {
            listData = data.results;
        } else if (data.data && Array.isArray(data.data)) {
            listData = data.data;
        } else {
            // 3. 리스트가 아예 없는 경우 (에러 메시지나 단일 객체일 수도 있음)
            console.warn("⚠️ 데이터가 리스트 형식이 아닙니다.");
            tbody.innerHTML = `<tr><td colspan="7" style="text-align:center; color:red;">
                Server Response is not a list. <br> Check Console for details.
            </td></tr>`;
            return;
        }
    }

    // 데이터가 비었으면 표시
    if (listData.length === 0) {
        tbody.innerHTML = '<tr><td colspan="7" style="text-align:center;">No candidates found.</td></tr>';
        return;
    }

    // 실제 그리기 (forEach는 이제 안전합니다)
    listData.forEach((item, idx) => {
        const tr = document.createElement("tr");
        tr.innerHTML = `
            <td>${idx + 1}</td>
            <td class="seq-display" style="font-size:11px;">${item.forward_primer || item.fwd_seq || '-'}</td>
            <td class="seq-display" style="font-size:11px;">${item.reverse_primer || item.rev_seq || '-'}</td>
            <td class="seq-display probe-color" style="font-size:11px;">${item.probe || item.probe_seq || '-'}</td>
            <td>${item.amplicon_size} bp</td>
            <td>${item.tm_f || '-'}/${item.tm_r || '-'}/${item.tm_p || '-'}</td>
            <td style="font-weight:bold;">${item.penalty ? item.penalty.toFixed(2) : '-'}</td>
        `;
        tbody.appendChild(tr);
    });
}
// 5. 메인 실행 로직
document.addEventListener("DOMContentLoaded", () => {
    
    // (A) 페이지 로드 시 저장된 값 복구
    restoreInputsFromStorage();

    const runBtn = document.getElementById("btn-run-qpcr");

    if (runBtn) {
        runBtn.addEventListener("click", async (e) => {
            e.preventDefault(); // 폼 전송에 의한 새로고침 방지

            // (B) 버튼 클릭 시 현재 값 저장
            saveInputsToStorage();

            // (C) 입력값 추출
            const payload = {
                // reference: "hg38", // 필요시 추가
                reference: document.getElementById("input-ref-genome").value,
                chrom: document.getElementById("input-chrom").value.trim(),
                start: parseInt(document.getElementById("input-start").value),
                end: parseInt(document.getElementById("input-end").value),
                ref: document.getElementById("input-ref-base")?.value.toUpperCase() || "G",
                alt: document.getElementById("input-alt-base")?.value.toUpperCase() || "A",
                top_k: 5,
                // 추가 옵션들...
                min_amplicon_length: parseInt(document.getElementById("min_amplicon_length")?.value) || 80,
                max_amplicon_length: parseInt(document.getElementById("max_amplicon_length")?.value) || 200
            };

            console.log("📤 [Web] Sending payload:", payload);

            // 유효성 검사
            if (!payload.chrom || isNaN(payload.start)) {
                alert("Please fill in Chromosome and Position.");
                return;
            }

            // (D) UI 로딩 상태 변경
            const originalBtnText = runBtn.innerText;
            runBtn.innerText = "⏳ RUNNING...";
            runBtn.disabled = true;
            if(document.getElementById("summary-status")) {
                document.getElementById("summary-status").innerText = "RUNNING";
            }

            try {
                // (E) 실제 FastAPI 서버 통신
                const response = await fetch("http://localhost:8080/api/design/qpcr", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify(payload)
                });

                if (!response.ok) {
                    const errData = await response.json();
                    throw new Error(errData.detail || "Server Error");
                }

                const results = await response.json();
                console.log("📥 [Web] Received results:", results);

                // (F) 결과 렌더링
                renderResults(results);

                // Summary 업데이트
                if(document.getElementById("summary-chrom")) document.getElementById("summary-chrom").innerText = payload.chrom;
                if(document.getElementById("summary-pos")) document.getElementById("summary-pos").innerText = `${payload.start}-${payload.end}`;
                if(document.getElementById("summary-status")) document.getElementById("summary-status").innerText = "COMPLETED";

            } catch (error) {
                console.error("Fetch Error:", error);
                alert("Design failed: " + error.message);
                if(document.getElementById("summary-status")) document.getElementById("summary-status").innerText = "FAILED";
            } finally {
                // (G) 버튼 원상복구
                runBtn.innerText = originalBtnText;
                runBtn.disabled = false;
            }
        });
    } else {
        console.warn("⚠️ 'btn-run-qpcr' 버튼을 찾을 수 없습니다. HTML ID를 확인하세요.");
    }
});