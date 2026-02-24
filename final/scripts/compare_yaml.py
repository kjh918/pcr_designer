#!/usr/bin/env python3
import sys
import os
import yaml
import argparse
from typing import Any, Dict

# 색상 코드 (가독성용)
GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
RESET = "\033[0m"
BLUE = "\033[94m"

def load_yaml(path: str) -> Dict[str, Any]:
    if not os.path.exists(path):
        print(f"{RED}❌ Error: File not found: {path}{RESET}")
        sys.exit(1)
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    except Exception as e:
        print(f"{RED}❌ Error parsing YAML: {e}{RESET}")
        sys.exit(1)

def compare_recursive(d1: Dict, d2: Dict, path: str = ""):
    """
    두 딕셔너리를 재귀적으로 비교합니다.
    d1: 기준 (Base)
    d2: 비교 대상 (Target/Temp)
    """
    all_keys = set(d1.keys()) | set(d2.keys())
    
    for key in sorted(all_keys):
        current_path = f"{path}.{key}" if path else key
        
        # 1. 한쪽에만 키가 있는 경우
        if key not in d1:
            print(f"{GREEN}[ADDED]{RESET}   {current_path}: {d2[key]}")
            continue
        if key not in d2:
            print(f"{RED}[MISSING]{RESET} {current_path} (was: {d1[key]})")
            continue
            
        val1 = d1[key]
        val2 = d2[key]

        # 2. 둘 다 딕셔너리인 경우 -> 재귀 호출
        if isinstance(val1, dict) and isinstance(val2, dict):
            compare_recursive(val1, val2, current_path)
        
        # 3. 값이 다른 경우
        elif val1 != val2:
            print(f"{YELLOW}[DIFF]{RESET}    {current_path}")
            print(f"   └─ Base: {val1}")
            print(f"   └─ Temp: {val2}")

def main():
    #parser = argparse.ArgumentParser(description="Compare two YAML files recursively.")
    #parser.add_argument("base_yaml", help="Path to the source/base YAML file")
    #parser.add_argument("/tmp/tmpxt9olbj5.yaml", help="Path to the target/temp YAML file")
    
    #args = parser.parse_args()

    #print(f"\n🔹 {BLUE}Comparing YAML Files{RESET}")
    #print(f"   Base: {args.base_yaml}")
    #print(f"   Temp: {args.temp_yaml}")
    #print("-" * 60)

    data1 = load_yaml('/home/jhkim/project/pcr_designer/final/pcr/config/base_pcr.yaml')
    data2 = load_yaml('/tmp/tmpxt9olbj5.yaml')

    compare_recursive(data1, data2)
    print("-" * 60)
    print("Done.")

if __name__ == "__main__":
    main()