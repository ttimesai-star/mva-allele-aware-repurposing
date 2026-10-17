"""Build the Track 2 report from report_src.txt.

- resolves [@Key; @Key2] citations to numbered references (order of first use) from config/references.yaml
  and outputs/reference_verification.csv (PubMed-verified metadata); an unknown key stops the build
- writes <name>.md (figures linked to the GitHub repository) and <name>.html (figures embedded from ../figures)
- optional: --pdf renders the HTML to PDF with a headless Chromium browser (Edge or Chrome)

python report/build_report.py --name zbignevich1_track2_report --pdf
"""
from __future__ import annotations

import argparse
import base64
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import markdown
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
GITHUB_FIG = "https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/"
CITE = re.compile(r"\[(@[A-Za-z0-9]+(?:\s*;\s*@[A-Za-z0-9]+)*)\]")

CSS = """
@page { size: A4; margin: 16mm 15mm 16mm 15mm; }
body { font-family: 'Segoe UI', Arial, sans-serif; font-size: 10pt; line-height: 1.42; color: #1f2933; }
h1 { font-size: 17pt; margin: 0 0 2pt 0; color: #102a43; }
.subtitle { font-size: 10.5pt; color: #52606d; margin-bottom: 10pt; }
h2 { font-size: 13pt; color: #102a43; border-bottom: 1px solid #bcccdc; padding-bottom: 2pt; margin-top: 16pt; }
h3 { font-size: 11pt; color: #243b53; margin-top: 12pt; }
blockquote { background: #fff8e6; border-left: 4px solid #b7791f; margin: 8pt 0; padding: 6pt 10pt; }
code { font-family: Consolas, monospace; font-size: 8.8pt; background: #f0f4f8; padding: 0 2px; }
table { border-collapse: collapse; font-size: 8.8pt; margin: 6pt 0; width: 100%; }
th, td { border: 1px solid #bcccdc; padding: 3pt 5pt; text-align: center; }
th { background: #e6edf4; } td:first-child, th:first-child { text-align: left; }
img { max-width: 100%; display: block; margin: 8pt auto 2pt auto; }
p.cap { font-size: 8.8pt; color: #52606d; margin-top: 0; }
ol.refs { font-size: 8.6pt; padding-left: 18pt; } ol.refs li { margin-bottom: 2pt; }
"""


def load_refs() -> dict:
    refs = yaml.safe_load((ROOT / "config" / "references.yaml").read_text(encoding="utf-8"))
    meta = {r["key"]: r for r in csv.DictReader(open(ROOT / "outputs" / "reference_verification.csv", encoding="utf-8"))}
    for k, v in refs.items():
        m = meta.get(k, {})
        if m.get("status") != "VERIFIED":
            sys.exit(f"reference {k} is not VERIFIED in outputs/reference_verification.csv")
        v.update({"journal": m.get("journal", ""), "year": m.get("year", ""), "author": m.get("first_author", ""),
                  "ptitle": m.get("pubmed_title") or v["title"]})
    return refs


def ref_line(r: dict) -> str:
    doi = f" doi:{r['doi']}." if r.get("doi") else ""
    return f"{r['author']} et al. {r['ptitle'].rstrip('.')}. *{r['journal']}* {r['year']}. PMID {r['pmid']}.{doi}"


def build(src: str, refs: dict) -> tuple[str, list[str]]:
    order: list[str] = []

    def sub(m):
        keys = [k.strip().lstrip("@") for k in m.group(1).split(";")]
        nums = []
        for k in keys:
            if k not in refs:
                sys.exit(f"unknown citation key: {k}")
            if k not in order:
                order.append(k)
            nums.append(order.index(k) + 1)
        return "[" + ", ".join(str(n) for n in sorted(set(nums))) + "]"

    body = CITE.sub(sub, src)
    lines = [f"{i}. {ref_line(refs[k])}" for i, k in enumerate(order, 1)]
    rest = [k for k in refs if k not in order]
    txt = "\n".join(lines)
    if rest:
        txt += ("\n\n*Additional references verified by the pipeline and used in `config/evidence.yaml`:*\n\n" +
                "\n".join(f"{len(order) + i}. {ref_line(refs[k])}" for i, k in enumerate(rest, 1)))
    return body.replace("{{REFERENCES}}", txt), order


def split_front(src: str) -> tuple[dict, str]:
    m = re.match(r"---\n(.*?)\n---\n", src, re.S)
    return (yaml.safe_load(m.group(1)), src[m.end():]) if m else ({}, src)


def to_html(front: dict, body: str) -> str:
    def img(m):
        cap, path = m.group(1), m.group(2).replace("FIG/", "")
        data = base64.b64encode((ROOT / "figures" / path).read_bytes()).decode()
        return f'<img src="data:image/png;base64,{data}" alt="figure"><p class="cap">{cap}</p>'

    body = re.sub(r"!\[(.*?)\]\((FIG/[^)]+)\)", img, body)
    html = markdown.markdown(body, extensions=["tables", "sane_lists"])
    html = html.replace("<ol>\n<li>" , "<ol>\n<li>", 1)
    html = re.sub(r'(<h2>References</h2>\s*)<ol>', r'\1<ol class="refs">', html)
    html = html.replace("<p><em>Additional references", '<p class="cap"><em>Additional references')
    html = re.sub(r'(Additional references.*?</p>\s*)<ol start="(\d+)">', r'\1<ol class="refs" start="\2">', html, flags=re.S)
    return (f"<!doctype html><html><head><meta charset='utf-8'><title>{front.get('title','')}</title>"
            f"<style>{CSS}</style></head><body><h1>{front.get('title','')}</h1>"
            f"<div class='subtitle'>{front.get('subtitle','')}</div>{html}</body></html>")


def browser() -> str | None:
    for p in [r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
              r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
              shutil.which("chromium") or "", shutil.which("google-chrome") or ""]:
        if p and os.path.exists(p):
            return p
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="zbignevich1_track2_report")
    ap.add_argument("--pdf", action="store_true")
    a = ap.parse_args()
    front, src = split_front((HERE / "report_src.txt").read_text(encoding="utf-8"))
    body, order = build(src, load_refs())
    md = f"# {front['title']}\n\n*{front['subtitle']}*\n\n" + re.sub(r"\]\(FIG/", f"]({GITHUB_FIG}", body)
    (HERE / f"{a.name}.md").write_text(md, encoding="utf-8")
    html_path = HERE / f"{a.name}.html"
    html_path.write_text(to_html(front, body), encoding="utf-8")
    print(f"cited {len(order)} references; wrote {a.name}.md and .html")
    if a.pdf:
        b = browser()
        if not b:
            sys.exit("no Chromium browser found for PDF rendering")
        with tempfile.TemporaryDirectory() as prof:  # throw-away profile: never touches the user's browser data
            subprocess.run([b, "--headless=new", "--disable-gpu", f"--user-data-dir={prof}", "--no-pdf-header-footer",
                            f"--print-to-pdf={HERE / (a.name + '.pdf')}", html_path.as_uri()], check=True, timeout=180)
        print(f"wrote {a.name}.pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
