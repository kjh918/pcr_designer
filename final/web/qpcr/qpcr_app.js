document.addEventListener("DOMContentLoaded", () => {
    const runBtn = document.getElementById("btn-run-qpcr");

    if (runBtn) {
        runBtn.addEventListener("click", async () => {
            // 1. Web UI에서 입력값 추출
            const payload = {
                chrom: document.getElementById("input-chrom").value.trim(),
                start: parseInt(document.getElementById("input-start").value),
                end: parseInt(document.getElementById("input-end").value),
                ref: document.getElementById("input-ref-base")?.value.toUpperCase() || "G",
                alt: document.getElementById("input-alt-base")?.value.toUpperCase() || "A",
                top_k: 5
            };

            // 유효성 검사
            if (!payload.chrom || isNaN(payload.start)) {
                alert("Please fill in Chromosome and Position.");
                return;
            }

            // 2. 버튼 로딩 상태 표시
            runBtn.innerText = "⏳ RUNNING...";
            runBtn.disabled = true;
            document.getElementById("summary-status").innerText = "RUNNING";

            try {
                // 3. FastAPI 서버로 요청
                const response = await fetch("http://localhost:8080/api/design/qpcr", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify(payload)
                });

                if (!response.ok) throw new Error("API Server Error");

                const results = await response.json();

                // 4. 결과 렌더링
                renderResults(results);
                
                // Summary 업데이트
                document.getElementById("summary-chrom").innerText = payload.chrom;
                document.getElementById("summary-pos").innerText = `${payload.start}-${payload.end}`;
                document.getElementById("summary-status").innerText = "COMPLETED";

            } catch (error) {
                console.error(error);
                alert("Design failed: " + error.message);
                document.getElementById("summary-status").innerText = "FAILED";
            } finally {
                runBtn.innerText = "DESIGN START";
                runBtn.disabled = false;
            }
        });
    }
});

/**
 * 서버에서 받은 JSON 데이터를 결과 테이블에 매핑합니다.
 */
function renderResults(data) {
    const tbody = document.getElementById("candidate-tbody");
    tbody.innerHTML = ""; 

    if (!data || data.length === 0) {
        tbody.innerHTML = '<tr><td colspan="7" style="text-align:center;">No candidates found.</td></tr>';
        return;
    }

    data.forEach((item, idx) => {
        const tr = document.createElement("tr");
        // 스크립트에서 리턴하는 JSON 키값에 맞춰 매핑
        tr.innerHTML = `
            <td>${idx + 1}</td>
            <td class="seq-display" style="font-size:11px;">${item.forward_primer || item.fwd_seq}</td>
            <td class="seq-display" style="font-size:11px;">${item.reverse_primer || item.rev_seq}</td>
            <td class="seq-display probe-color" style="font-size:11px;">${item.probe || item.probe_seq}</td>
            <td>${item.amplicon_size} bp</td>
            <td>${item.tm_f || '-'}/${item.tm_r || '-'}/${item.tm_p || '-'}</td>
            <td style="font-weight:bold;">${item.penalty ? item.penalty.toFixed(2) : '-'}</td>
        `;
        tbody.appendChild(tr);
    });
}

// 💡 추적할 입력창들의 ID 목록입니다. (HTML에 적어둔 id와 동일해야 함)
const TRACKED_INPUTS = [
    "input-ref-genome",
    "input-chrom",
    "input-start",
    "input-end",
    "input-strand",
    // 만약 qpcr_params.html에 아래 ID들이 있다면 같이 추적합니다.
    "input-ref-base", 
    "input-alt-base",
    "input-prod-size"
];

// ---------------------------------------------------------
// 1. 브라우저 창고(localStorage)에 값 저장하기
// ---------------------------------------------------------
function saveInputsToStorage() {
    TRACKED_INPUTS.forEach(id => {
        const inputEl = document.getElementById(id);
        if (inputEl && inputEl.value) {
            // 브라우저 창고에 'ID : 값' 형태로 저장
            localStorage.setItem("qpcr_" + id, inputEl.value); 
        }
    });
    console.log("✅ 현재 입력값이 브라우저에 안전하게 저장되었습니다!");
}

// ---------------------------------------------------------
// 2. 브라우저 창고에서 값 꺼내와서 채워넣기
// ---------------------------------------------------------
function restoreInputsFromStorage() {
    TRACKED_INPUTS.forEach(id => {
        const savedValue = localStorage.getItem("qpcr_" + id);
        if (savedValue !== null) {
            const inputEl = document.getElementById(id);
            if (inputEl) {
                inputEl.value = savedValue; // 입력창에 값 복구
            }
        }
    });
}

// ---------------------------------------------------------
// 3. 메인 실행 로직 (페이지가 켜지면 시작됨)
// ---------------------------------------------------------
document.addEventListener("DOMContentLoaded", () => {
    
    // 페이지가 로드되자마자 창고에서 예전 기록을 꺼내서 폼에 채워줍니다.
    restoreInputsFromStorage();

    const runBtn = document.getElementById("btn-run-qpcr");
    
    // [DESIGN START] 버튼을 클릭했을 때의 동작
    if (runBtn) {
        runBtn.addEventListener("click", () => {
            
            // 1) 클릭하는 순간 현재 적혀있는 값들을 영원히(?) 기억하게 저장
            saveInputsToStorage();

            // 2) UI를 "로딩 중" 상태로 변경
            runBtn.innerText = "⏳ RUNNING...";
            runBtn.style.backgroundColor = "#6c757d"; // 회색으로 변경
            runBtn.disabled = true;

            // 3) 파이썬 백엔드(FastAPI)로 데이터 전송 (현재는 2초 대기하는 가짜 로직)
            setTimeout(() => {
                alert("디자인이 완료되었습니다! (임시 메시지)");
                
                // 버튼 원상복구
                runBtn.innerText = "DESIGN START";
                runBtn.style.backgroundColor = "#0056b3";
                runBtn.disabled = false;
                
                // TODO: 여기서 백엔드에서 받은 결과를 테이블(result_panel)에 그려줍니다.
            }, 2000);
            
        });
    }
});