"""Allele-class landscape of reported pathogenic variants in MVA genes (ClinVar, public).

For each gene: pull Pathogenic / Likely pathogenic single-gene variants, parse the protein
change from the ClinVar title, assign an allele class (config/evidence.yaml) and report how
many variants fall in each class. For BUB1B this tells us what fraction of MVA1 alleles
could respond to a dosage-raising strategy (S1) at all.

No hackathon (controlled) data is used anywhere in this repository.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import Counter

from common import OUTPUTS, PROCESSED, eutils, fetch_json, load_yaml

GENES = ["BUB1B", "TRIP13", "CEP57", "BUB1", "MAD1L1", "CENATAC", "BUB3"]
UNIPROT = {"BUB1B": "O60566"}
AA3 = "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val".split()


def classify(title: str, consequences: list[str], window: tuple[int, int]) -> tuple[str, int | None]:
    """Return (allele_class, protein_position)."""
    p = re.search(r"\(p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*|fs|=)?([^)]*)\)", title)
    pos = int(p.group(2)) if p else None
    t = title.lower()
    cons = " ".join(consequences).lower()
    if p and (p.group(3) in ("Ter", "*") or "fs" in (p.group(3) or "") + p.group(4) or "ter" in p.group(4).lower()):
        return "TRUNC", pos
    if re.search(r"c\.\d+[+-][12][acgt]>", t) or "splice donor" in cons or "splice acceptor" in cons:
        return "TRUNC", pos
    if "frameshift" in cons or "nonsense" in cons or re.search(r"c\.[^ ]*(del|dup|ins)", t) and not p:
        return "TRUNC", pos
    if p and p.group(3) in AA3 and p.group(1) in AA3:
        return ("MS_UNSTABLE" if window[0] <= pos <= window[1] else "MS_FUNCTIONAL"), pos
    if re.search(r"c\.-\d+|upstream|5'utr|intergenic", t + " " + cons):
        return "REG", pos
    return "OTHER", pos


def gene_variants(gene: str) -> list[dict]:
    term = f'{gene}[gene] AND (clinsig_pathogenic[prop] OR clinsig_likely_pathogenic[prop]) AND single_gene[prop]'
    ids = eutils("esearch", db="clinvar", term=term, retmax=1000)["esearchresult"]["idlist"]
    out = []
    for i in range(0, len(ids), 100):
        res = eutils("esummary", db="clinvar", id=",".join(ids[i:i + 100]))["result"]
        for uid in res.get("uids", []):
            v = res[uid]
            title = v.get("title", "")
            if f"({gene})" not in title:
                continue
            gc = v.get("germline_classification", {}) or {}
            traits = "; ".join(t.get("trait_name", "") for t in gc.get("trait_set", []) or [])
            out.append({"gene": gene, "clinvar_uid": uid, "title": title,
                        "classification": gc.get("description", ""), "review": gc.get("review_status", ""),
                        "traits": traits, "consequences": v.get("molecular_consequence_list", []) or []})
    return out


def main() -> int:
    ev = load_yaml("evidence.yaml")
    win = (ev["kinase_window"]["start"], ev["kinase_window"]["end"])
    rows, summary = [], []
    for g in GENES:
        vs = gene_variants(g)
        for v in vs:
            cls, pos = classify(v["title"], v["consequences"], win if g == "BUB1B" else (10**9, 10**9))
            if cls == "TRUNC" and re.search(r"Ter\)|\*\)", v["title"]) and "fs" not in v["title"]:
                cls = "TRUNC_NONSENSE"
            if g != "BUB1B" and cls == "MS_FUNCTIONAL":
                cls = "MISSENSE"
            v.update({"allele_class": cls, "protein_pos": pos, "consequences": "|".join(v["consequences"])})
            rows.append(v)
        c = Counter(r["allele_class"] for r in rows if r["gene"] == g)
        summary.append({"gene": g, "n_pathogenic_or_LP": len(vs), **dict(c)})
        print(g, len(vs), dict(c))
    fields = ["gene", "clinvar_uid", "title", "classification", "review", "traits", "consequences", "allele_class", "protein_pos"]
    with open(PROCESSED / "clinvar_mva_variants.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    keys = sorted({k for s in summary for k in s} - {"gene", "n_pathogenic_or_LP"})
    with open(OUTPUTS / "clinvar_allele_class_summary.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "n_pathogenic_or_LP"] + keys)
        w.writeheader()
        for s in summary:
            w.writerow({k: s.get(k, 0) for k in ["gene", "n_pathogenic_or_LP"] + keys})
    # How are BUB1B missense variants with an MVA/PCS trait classified in ClinVar?
    term = ('BUB1B[gene] AND missense_variant[molcons] AND single_gene[prop] AND '
            '("mosaic variegated aneuploidy"[dis] OR "premature chromatid separation"[dis])')
    ids = eutils("esearch", db="clinvar", term=term, retmax=1000)["esearchresult"]["idlist"]
    cl = Counter()
    for i in range(0, len(ids), 200):
        res = eutils("esummary", db="clinvar", id=",".join(ids[i:i + 200]))["result"]
        for uid in res.get("uids", []):
            cl[(res[uid].get("germline_classification") or {}).get("description", "")] += 1
    with open(OUTPUTS / "bub1b_missense_clinvar_classes.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["classification", "n"])
        for k, n in cl.most_common():
            w.writerow([k, n])
    print("BUB1B missense with MVA/PCS trait:", dict(cl))
    # Literature check of the kinase-window rule (Suijkerbuijk 2010, PMC2887387)
    lit = list(csv.DictReader(open(PROCESSED.parent / "literature_mva1_missense.csv", encoding="utf-8")))
    agree = 0
    for r in lit:
        pos = int(re.search(r"\d+", r["variant"]).group())
        pred = "MS_UNSTABLE" if win[0] <= pos <= win[1] else "MS_FUNCTIONAL"
        r["predicted_class"] = pred
        r["agrees"] = (pred == "MS_UNSTABLE") == (r["observed_abundance"] == "low")
        agree += r["agrees"]
    with open(OUTPUTS / "kinase_window_rule_check.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(lit[0]))
        w.writeheader()
        w.writerows(lit)
    print(f"kinase-window rule vs published abundance data: {agree}/{len(lit)} concordant")
    # UniProt domain map for the lollipop figure
    up = fetch_json(f"https://rest.uniprot.org/uniprotkb/{UNIPROT['BUB1B']}.json?fields=ft_domain,ft_region,ft_motif,length")
    with open(PROCESSED / "bubr1_domains.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["type", "start", "end", "description"])
        for f in up.get("features", []):
            w.writerow([f["type"], f["location"]["start"]["value"], f["location"]["end"]["value"], f.get("description", "")])
    return 0


if __name__ == "__main__":
    sys.exit(main())
