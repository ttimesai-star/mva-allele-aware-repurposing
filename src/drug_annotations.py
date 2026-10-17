"""Live annotation of candidate drugs + automatic do-no-harm scan (ChEMBL, openFDA, Open Targets).

Outputs
  outputs/candidate_annotations.csv  - ChEMBL id, max phase, first approval, mechanisms, pediatric label text
  outputs/liability_scan.csv         - every approved (max_phase 4) drug whose ChEMBL mechanism matches a
                                       do-no-harm pattern in evidence.yaml
  outputs/mva_gene_associations.csv  - Open Targets genes associated with MVA (scalability input)
"""
from __future__ import annotations

import csv
import sys
import urllib.parse

from common import OUTPUTS, fetch_json, load_yaml

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"


def chembl_molecule(name: str) -> dict:
    url = f"{CHEMBL}/molecule.json?pref_name__iexact={urllib.parse.quote(name)}&limit=5"
    mols = fetch_json(url).get("molecules", [])
    if not mols:  # fall back to synonym search
        url = f"{CHEMBL}/molecule/search.json?q={urllib.parse.quote(name)}&limit=5"
        mols = fetch_json(url).get("molecules", [])
    mols.sort(key=lambda m: -(float(m.get("max_phase") or 0)))
    return mols[0] if mols else {}


def chembl_mechanisms(chembl_id: str) -> list[str]:
    url = f"{CHEMBL}/mechanism.json?molecule_chembl_id={chembl_id}&limit=50"
    return sorted({m.get("mechanism_of_action", "") for m in fetch_json(url).get("mechanisms", [])})


def openfda_pediatric(name: str) -> str:
    q = urllib.parse.quote(f'openfda.generic_name:"{name}"')
    try:
        res = fetch_json(f"https://api.fda.gov/drug/label.json?search={q}&limit=25").get("results", [])
    except RuntimeError:
        return ""
    texts = []
    for r in res:
        gn = [g.lower() for g in r.get("openfda", {}).get("generic_name", [])]
        if gn and all(name.lower() in g and " and " not in g for g in gn) and r.get("pediatric_use"):
            texts.append(r["pediatric_use"][0])
    texts.sort(key=len, reverse=True)
    return " ".join(texts[0].split())[:600] if texts else ""


def liability_scan(patterns: list[dict]) -> list[dict]:
    rows = []
    for p in patterns:
        url = f"{CHEMBL}/mechanism.json?mechanism_of_action__icontains={urllib.parse.quote(p['pattern'])}&limit=1000"
        mechs = fetch_json(url).get("mechanisms", [])
        ids = sorted({m["molecule_chembl_id"] for m in mechs})
        moa = {m["molecule_chembl_id"]: m.get("mechanism_of_action", "") for m in mechs}
        for i in range(0, len(ids), 50):
            chunk = ",".join(ids[i:i + 50])
            mols = fetch_json(f"{CHEMBL}/molecule.json?molecule_chembl_id__in={chunk}&limit=50").get("molecules", [])
            for m in mols:
                if float(m.get("max_phase") or 0) >= 4:
                    rows.append({"pattern": p["pattern"], "drug": m.get("pref_name"), "chembl_id": m["molecule_chembl_id"],
                                 "mechanism": moa.get(m["molecule_chembl_id"], ""), "first_approval": m.get("first_approval"),
                                 "why_flagged": p["why"], "refs": ";".join(p["refs"])})
    # one row per drug (a drug can match two patterns)
    seen, out = set(), []
    for r in sorted(rows, key=lambda r: (r["pattern"], r["drug"] or "")):
        if r["drug"] and r["drug"] not in seen:
            seen.add(r["drug"])
            out.append(r)
    return out


def ot_mva_genes() -> list[dict]:
    q = """query($id:String!){disease(efoId:$id){name associatedTargets(page:{index:0,size:40}){rows{score target{approvedSymbol approvedName}}}}}"""
    d = fetch_json("https://api.platform.opentargets.org/api/v4/graphql", body={"query": q, "variables": {"id": "MONDO_0000141"}})
    rows = d["data"]["disease"]["associatedTargets"]["rows"]
    return [{"gene": r["target"]["approvedSymbol"], "name": r["target"]["approvedName"], "ot_score": round(r["score"], 3)}
            for r in rows if r["score"] >= 0.3]


def main() -> int:
    ev = load_yaml("evidence.yaml")
    out, problems = [], []
    for c in ev["candidates"]:
        m = chembl_molecule(c["drug"])
        cid = m.get("molecule_chembl_id", "")
        parent = (m.get("molecule_hierarchy") or {}).get("parent_chembl_id", "")
        if c["chembl"] not in (cid, parent):
            problems.append(f"{c['drug']}: evidence.yaml says {c['chembl']}, ChEMBL resolves {cid}")
        out.append({"drug": c["drug"], "chembl_id": cid, "pref_name": m.get("pref_name"), "max_phase": m.get("max_phase"),
                    "first_approval": m.get("first_approval"), "mechanisms": "; ".join(chembl_mechanisms(cid)) if cid else "",
                    "fda_pediatric_use": openfda_pediatric(c["drug"])})
        print(f"{c['drug']:24s} {cid:14s} phase={m.get('max_phase')} first={m.get('first_approval')}")
    with open(OUTPUTS / "candidate_annotations.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    lia = liability_scan(ev["liability_mechanisms"])
    with open(OUTPUTS / "liability_scan.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(lia[0]))
        w.writeheader()
        w.writerows(lia)
    print(f"do-no-harm scan: {len(lia)} approved drugs flagged")
    genes = ot_mva_genes()
    with open(OUTPUTS / "mva_gene_associations.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(genes[0]))
        w.writeheader()
        w.writerows(genes)
    print("Open Targets MVA genes (score>=0.3):", ", ".join(g["gene"] for g in genes))
    for p in problems:
        print("CHEMBL ID MISMATCH:", p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
