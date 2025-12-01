# mypkg/core.py

from typing import Dict

def run_analysis(sample_name: str, threshold: float) -> Dict:
    # TODO: 여기에 실제 패키지 로직을 넣으면 됩니다.
    # 아래는 그냥 형태 예시
    score = len(sample_name) * threshold
    return {
        "sample_name": sample_name,
        "threshold": threshold,
        "score": score,
        "message": f"{sample_name} 분석 완료"
    }
