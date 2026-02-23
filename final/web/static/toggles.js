<script>
document.addEventListener("DOMContentLoaded", () => {
  // 버튼 이벤트 리스너 등록
  const runBtn = document.querySelector(".btn-run");
  if(runBtn) {
    runBtn.addEventListener("click", runDesignPipeline);
  }
});

async function runDesignPipeline() {
  const statusEl = document.getElementById("summary-status");
  statusEl.textContent = "RUNNING...";
  statusEl.style.color = "blue";
  
  try {
    // 💡 실제 환경에서는 FastAPI/Flask 등 백엔드 API 주소로 fetch 해야 합니다.
    // 예: const response = await fetch("http://localhost:8000/api/design", { method: "POST", body: ... });
    // 여기서는 통신이 성공하여 방금 파이썬이 만든 JSON 객체(data)를 받았다고 가정합니다.
    
    // 임시 테스트용 목업 데이터 (파이썬 스크립트의 결과물 구조와 동일)
    const data = await getMockData(); 

    if (data.status !== "success") {
      statusEl.textContent = "FAILED: " + data.reason;
      statusEl.style.color = "red";
      return;
    }

    // 1. 요약(Summary) 업데이트
    statusEl.textContent = "SUCCESS (" + data.metadata.candidates_found + " found)";
    statusEl.style.color = "green";
    
    const topResult = data.results[0];
    document.getElementById("summary-size").textContent = topResult.amplicon_info.length + " bp";

    // 2. Alignment View 업데이트
    document.getElementById("alignment-view").textContent = topResult.alignment_text_block;

    // 3. 테이블(Candidate List) 렌더링
    const tbody = document.getElementById("candidate-tbody");
    tbody.innerHTML = ""; // 기존 내용 비우기

    data.results.forEach((res, index) => {
      const tr = document.createElement("tr");
      // 행을 클릭하면 하단 상세 뷰가 바뀌도록 이벤트 추가
      tr.style.cursor = "pointer";
      tr.addEventListener("click", () => updateDetailView(res));

      tr.innerHTML = `
        <td>${res.rank}</td>
        <td>${res.oligos.forward.sequence}</td>
        <td>${res.oligos.reverse.sequence}</td>
        <td>${res.oligos.probe ? res.oligos.probe.sequence : '-'}</td>
        <td>${res.amplicon_info.length} bp</td>
        <td>${res.oligos.forward.tm} / ${res.oligos.reverse.tm} / ${res.oligos.probe ? res.oligos.probe.tm : '-'}</td>
        <td>${res.oligos.forward.gc} / ${res.oligos.reverse.gc} / ${res.oligos.probe ? res.oligos.probe.gc : '-'}</td>
        <td>${res.metrics.pair_penalty}</td>
      `;
      tbody.appendChild(tr);
    });

    // 4. 초기 상세 뷰 세팅 (Rank 1 기준)
    updateDetailView(topResult);

  } catch (error) {
    console.error(error);
    statusEl.textContent = "ERROR";
    statusEl.style.color = "red";
  }
}

// 하단 상세 패널 및 QC 플래그 업데이트 함수
function updateDetailView(result) {
  // 서열 및 좌표 매핑
  document.getElementById("detail-fwd-seq").textContent = result.oligos.forward.sequence;
  document.getElementById("detail-fwd-pos").textContent = result.oligos.forward.genomic_pos;
  
  document.getElementById("detail-rev-seq").textContent = result.oligos.reverse.sequence;
  document.getElementById("detail-rev-pos").textContent = result.oligos.reverse.genomic_pos;
  
  if(result.oligos.probe) {
    document.getElementById("detail-prb-seq").textContent = result.oligos.probe.sequence;
    document.getElementById("detail-prb-pos").textContent = result.oligos.probe.genomic_pos;
  } else {
    document.getElementById("detail-prb-seq").textContent = "N/A";
    document.getElementById("detail-prb-pos").textContent = "";
  }

  // QC 로직 매핑
  const qcList = document.getElementById("qc-flags-list");
  qcList.innerHTML = "";
  
  const isPass = result.qc_info.is_pass;
  const qcColor = isPass ? "green" : "red";
  const qcIcon = isPass ? "✅" : "❌";
  const reasonText = isPass ? "All Passed" : result.qc_info.fail_reason;

  qcList.innerHTML = `
    <li><span style="color:${qcColor}; font-weight:bold;">Overall QC: ${qcIcon} ${reasonText}</span></li>
    <li>Amplicon Position: ${result.amplicon_info.genomic_pos}</li>
    <li>Amplicon Tm: ${result.amplicon_info.tm} °C</li>
    <li>Amplicon GC: ${result.amplicon_info.gc} %</li>
  `;
  
  // Alignment View 연동 변경 (클릭한 행에 맞춰 다이어그램도 바뀜)
  document.getElementById("alignment-view").textContent = result.alignment_text_block;
}

// ---------------------------------------------------
// 테스트용 가짜 데이터 제너레이터 (백엔드 연결 전 UI 확인용)
// ---------------------------------------------------
function getMockData() {
  return Promise.resolve({
    "status": "success",
    "metadata": { "candidates_found": 1 },
    "results": [
      {
        "rank": 1,
        "metrics": { "pair_penalty": 3.114 },
        "qc_info": { "is_pass": true, "fail_reason": "None" },
        "amplicon_info": {
          "length": 83,
          "tm": 78.45,
          "gc": 57.8,
          "genomic_pos": "chr17:39723541-39723623"
        },
        "oligos": {
          "forward": { "sequence": "GCAGATGCGGATCCTG", "tm": 54.35, "gc": 62.5, "genomic_pos": "chr17:39723541-39723556" },
          "reverse": { "sequence": "TGACCTTGTAGACTGTGC", "tm": 53.8, "gc": 55.6, "genomic_pos": "chr17:39723606-39723623" },
          "probe": { "sequence": "AGACGGAGCTGAGGAAGGTGAAGGG", "tm": 61.43, "gc": 64.0, "genomic_pos": "chr17:39723570-39723594" }
        },
        "alignment_text_block": "AMPLICON : GCAGATGCGGATCCTG...AGACGGAGCTGAGGAAGGTGAAGGG...GCACAGTCTACAAGGTCA\nFORWARD  : GCAGATGCGGATCCTG\nPROBE(+) :                   AGACGGAGCTGAGGAAGGTGAAGGG\nREVERSE  :                                               GCACAGTCTACAAGGTCA (RC)"
      }
    ]
  });
}
</script>