"""Turn a BUB1B genotype (two HGVS protein or cDNA descriptions) into allele classes and the
matching ranked list of candidate medicines.

Usage (run locally by whoever holds the confirmed variants; nothing is sent anywhere):
  python src/classify_variants.py "c.1833delT (p.Thr612fs)" "p.Ile909Thr"
  python src/classify_variants.py --file my_variants.txt      # one variant per line, two lines

This is how the hackathon panel (or a treating team) can instantiate the framework with the
clinically confirmed genotype without that genotype ever entering this repository.
"""
from __future__ import annotations

import argparse
import re
import sys

from clinvar_landscape import classify
from common import load_yaml
from score import base_score, fit


def to_title(v: str) -> str:
    v = v.strip()
    if v.startswith("p.") and "(" not in v:
        v = f"({v})"
    return v


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("variants", nargs="*")
    ap.add_argument("--file")
    a = ap.parse_args(argv)
    vs = a.variants or [l for l in open(a.file, encoding="utf-8").read().splitlines() if l.strip()]
    if len(vs) != 2:
        print("need exactly two alleles (biallelic MVA1)")
        return 2
    ev, cfg = load_yaml("evidence.yaml"), load_yaml("weights.yaml")
    win = (ev["kinase_window"]["start"], ev["kinase_window"]["end"])
    classes = []
    for v in vs:
        title = to_title(v)
        p = re.search(r"\(p\.", title)
        cls, pos = classify(title if p else title, [], win)
        if cls == "TRUNC" and re.search(r"Ter|\*", title) and "fs" not in title:
            cls = "TRUNC_NONSENSE"
        if cls == "OTHER":
            print(f"  ! could not classify '{v}' - give the protein change, e.g. p.Ile909Thr")
        classes.append(cls)
        print(f"{v:40s} -> {cls} (protein position {pos})")
    print("\nRanked candidates for this genotype:")
    rows = []
    for c in ev["candidates"]:
        st = ev["strategies"][c["strategy"]]
        f = fit(st, classes)
        b = base_score(c["scores"], cfg["weights"])
        rows.append((b * (cfg["fit_floor"] + (1 - cfg["fit_floor"]) * f), c["drug"], c["strategy"], f))
    for i, (s, d, st, f) in enumerate(sorted(rows, reverse=True), 1):
        print(f"{i:2d}. {d:24s} {st:20s} fit={f:.1f} score={s:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
