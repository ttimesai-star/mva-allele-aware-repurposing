"""Verify every reference in config/references.yaml against PubMed.

A reference passes when (a) the PMID resolves, (b) >= 70% of the content words of the
stored title occur in the PubMed title, and (c) the stored DOI (if any) equals the PubMed DOI.
Also checks that every reference key cited in evidence.yaml exists here.
Exit code 1 if anything fails, so `make all` stops on a broken citation.
"""
from __future__ import annotations

import csv
import re
import sys
import unicodedata

from common import OUTPUTS, eutils, load_yaml

STOP = {"a", "an", "the", "of", "in", "and", "with", "to", "by", "for", "on", "from", "is", "at", "as", "that", "are"}


def words(s: str) -> set[str]:
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode().lower()
    return {w for w in re.findall(r"[a-z0-9]+", s) if w not in STOP and len(w) > 1}


def cited_keys(obj) -> set[str]:
    out: set[str] = set()
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in ("refs", "counter_refs") and isinstance(v, list):
                out.update(v)
            else:
                out |= cited_keys(v)
    elif isinstance(obj, list):
        for v in obj:
            out |= cited_keys(v)
    return out


def main() -> int:
    refs = load_yaml("references.yaml")
    ids = [r["pmid"] for r in refs.values()]
    summ = {}
    for i in range(0, len(ids), 50):
        res = eutils("esummary", db="pubmed", id=",".join(ids[i:i + 50]))["result"]
        summ.update({k: v for k, v in res.items() if k != "uids"})
    rows, bad = [], 0
    for key, r in refs.items():
        s = summ.get(r["pmid"], {})
        title = s.get("title", "")
        doi = next((a["value"] for a in s.get("articleids", []) if a["idtype"] == "doi"), "")
        want = words(r["title"])
        overlap = len(want & words(title)) / max(1, len(want))
        doi_ok = (not r["doi"]) or (r["doi"].lower() == doi.lower())
        ok = bool(title) and overlap >= 0.7 and doi_ok
        bad += not ok
        rows.append({"key": key, "pmid": r["pmid"], "doi": doi, "journal": s.get("source", ""),
                     "year": s.get("pubdate", "")[:4], "first_author": (s.get("authors") or [{"name": ""}])[0]["name"],
                     "pubmed_title": title, "title_overlap": round(overlap, 2), "doi_match": doi_ok,
                     "status": "VERIFIED" if ok else "FAILED"})
    ev = load_yaml("evidence.yaml")
    missing = sorted(cited_keys(ev) - set(refs))
    with open(OUTPUTS / "reference_verification.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    n_ok = sum(r["status"] == "VERIFIED" for r in rows)
    print(f"references verified: {n_ok}/{len(rows)}; cited-but-undefined keys: {missing or 'none'}")
    for r in rows:
        if r["status"] != "VERIFIED":
            print("  FAILED", r)
    return 1 if (bad or missing) else 0


if __name__ == "__main__":
    sys.exit(main())
