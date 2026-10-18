"""Render a Markdown document to a print-ready PDF.

Converts the Markdown to a styled standalone HTML, then drives a headless Chromium-based browser
(Microsoft Edge or Google Chrome, whichever is installed) to print it. This avoids a heavy PDF
toolchain (LaTeX/pandoc/wkhtmltopdf) that is not available on this machine.

Usage:
    python tools/md2pdf.py MUDANCAS.md MUDANCAS.pdf
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import markdown

BROWSERS = [
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
]

CSS = """
@page { size: A4; margin: 18mm 16mm; }
* { box-sizing: border-box; }
body {
  font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
  font-size: 10.5pt; line-height: 1.55; color: #1b1f24; margin: 0;
}
h1 {
  font-size: 21pt; margin: 0 0 4pt; color: #0f2c4d; letter-spacing: -.2pt;
  border-bottom: 3px solid #0f2c4d; padding-bottom: 8pt;
}
h2 {
  font-size: 14pt; margin: 22pt 0 7pt; color: #0f2c4d;
  border-bottom: 1px solid #d4dbe3; padding-bottom: 4pt;
  page-break-after: avoid; break-after: avoid;
}
h3 { font-size: 11.5pt; margin: 14pt 0 5pt; color: #24425f; page-break-after: avoid; break-after: avoid; }
p { margin: 6pt 0; text-align: justify; }
strong { color: #0f2c4d; }
hr { border: 0; border-top: 1px solid #e2e7ec; margin: 16pt 0; }
ul, ol { margin: 6pt 0; padding-left: 20pt; }
li { margin: 3pt 0; }
code {
  font-family: "Cascadia Mono", Consolas, monospace; font-size: 9pt;
  background: #f2f5f8; padding: 1pt 4pt; border-radius: 3px; color: #10405f;
}
pre {
  background: #f7f9fb; border: 1px solid #e2e7ec; border-left: 3px solid #0f2c4d;
  border-radius: 4px; padding: 9pt 11pt; overflow-x: auto;
  page-break-inside: avoid; break-inside: avoid;
}
pre code { background: none; padding: 0; font-size: 8.7pt; line-height: 1.45; color: #1b1f24; }
table {
  border-collapse: collapse; width: 100%; margin: 9pt 0; font-size: 9pt;
  page-break-inside: avoid; break-inside: avoid;
}
th {
  background: #0f2c4d; color: #fff; text-align: left; font-weight: 600;
  padding: 5pt 7pt; border: 1px solid #0f2c4d;
}
td { padding: 5pt 7pt; border: 1px solid #d4dbe3; vertical-align: top; }
tr:nth-child(even) td { background: #f7f9fb; }
blockquote {
  margin: 9pt 0; padding: 7pt 12pt; background: #fdf8ec;
  border-left: 3px solid #c8901f; color: #4a3a15;
  page-break-inside: avoid; break-inside: avoid;
}
blockquote p { margin: 3pt 0; }
a { color: #10527f; text-decoration: none; }
em { color: #4a5560; }
"""

HTML = """<!doctype html>
<html lang="pt-BR"><head><meta charset="utf-8"><title>{title}</title>
<style>{css}</style></head><body>{body}</body></html>"""


def find_browser() -> str:
    for path in BROWSERS:
        if os.path.exists(path):
            return path
    raise SystemExit("No Chromium-based browser found for PDF printing (looked for Edge/Chrome).")


def render(md_path: Path, pdf_path: Path) -> Path:
    body = markdown.markdown(
        md_path.read_text(encoding="utf-8"),
        extensions=["tables", "fenced_code", "sane_lists", "attr_list"],
    )
    html = HTML.format(title=md_path.stem, css=CSS, body=body)

    with tempfile.TemporaryDirectory() as tmp:
        html_path = Path(tmp) / "doc.html"
        html_path.write_text(html, encoding="utf-8")
        subprocess.run(
            [
                find_browser(), "--headless", "--disable-gpu", "--no-sandbox",
                "--no-pdf-header-footer", f"--print-to-pdf={pdf_path}",
                html_path.as_uri(),
            ],
            check=True, capture_output=True, timeout=180,
        )
    if not pdf_path.exists():
        raise SystemExit(f"Browser did not produce {pdf_path}")
    return pdf_path


if __name__ == "__main__":
    src = Path(sys.argv[1] if len(sys.argv) > 1 else "MUDANCAS.md").resolve()
    dst = Path(sys.argv[2] if len(sys.argv) > 2 else src.with_suffix(".pdf")).resolve()
    out = render(src, dst)
    print(f"{out}  ({out.stat().st_size / 1024:.0f} KB)")
