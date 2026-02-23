# web/build_web.py
from pathlib import Path
import shutil
import subprocess

WEB = Path(__file__).resolve().parent
DIST = WEB / "dist"
QMD = WEB / "main.qmd"

def main() -> None:
    DIST.mkdir(parents=True, exist_ok=True)

    # ✅ web 디렉토리에서 렌더 (partials include 상대경로 안전)
    subprocess.check_call(["quarto", "render", QMD.name], cwd=str(WEB))

    # quarto output: web/main.html, web/main_files/
    out_html = WEB / "main.html"
    out_files = WEB / "main_files"

    # ✅ dist/index.html로 고정
    shutil.copy2(out_html, DIST / "index.html")

    # ✅ 리소스 폴더도 dist로 이동/동기화
    dest_files = DIST / "main_files"
    if out_files.exists():
        if dest_files.exists():
            shutil.rmtree(dest_files)
        shutil.copytree(out_files, dest_files)

    print(f"Built: {(DIST / 'index.html').resolve()}")

if __name__ == "__main__":
    main()