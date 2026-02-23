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