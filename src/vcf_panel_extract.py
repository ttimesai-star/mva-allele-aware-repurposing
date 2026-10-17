"""Extract and annotate variants in spindle-assembly-checkpoint / MVA genes from a single-sample VCF.

Run this ONLY on the machine that holds the controlled VCF. It writes its intermediate files to --workdir
(outside this repository) and prints a short summary: the rare, protein-altering or splice-relevant variants
in the gene panel. Nothing genome-wide is written or printed.

    python src/vcf_panel_extract.py --vcf /secure/WGS.vcf.gz --workdir /secure/work

Steps
1. Gene panel (16 SAC/MVA genes) -> MANE Select transcript, exons and CDS from the UCSC Genome Browser API (GRCh38).
2. Stream the VCF once; keep records within exons +/-20 bp and 2 kb upstream of each panel gene
   (intronic records are only counted).
3. Annotate those records (and only those) with Ensembl VEP REST: MANE/canonical consequence, HGVS, gnomAD
   exome/genome AF, AlphaMissense, REVEL, CADD, SpliceAI. Per the hackathon's data-use guidance an annotation API
   that returns a result and takes no rights in the data is a processing tool; we still send only panel exon records
   (about a hundred positions), never the genome.
4. Report rare (gnomAD max AF < 1% or absent) variants that are HIGH/MODERATE impact or splice-relevant
   (SpliceAI delta > 0.2).

The VCF used in the hackathon has contigs without the "chr" prefix and singleton genotypes; PGT/PID give
read-backed phase when two variants share a read pair.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import ssl
import time
import urllib.parse
import urllib.request

PANEL = ["BUB1B", "CEP57", "TRIP13", "CENATAC", "BUB1", "BUB3", "MAD1L1", "MAD2L1", "TTK", "CDC20", "KNL1",
         "ZWINT", "CENPE", "SGO1", "KIF2C", "ANAPC1"]
VEP = ("https://rest.ensembl.org/vep/homo_sapiens/region?canonical=1&hgvs=1&mane=1&domains=1&AlphaMissense=1"
       "&CADD=1&SpliceAI=1&REVEL=1&protein=1&af_gnomadg=1&af_gnomade=1&variant_class=1")


try:  # some Windows Pythons lack a CA bundle for urllib; certifi (shipped with requests) fixes it
    import certifi
    _CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    _CTX = ssl.create_default_context()


def _get(url: str, data: bytes | None = None, tries: int = 8):
    hdr = {"Content-Type": "application/json", "Accept": "application/json"}
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=hdr), timeout=180, context=_CTX) as r:
                return json.loads(r.read().decode())
        except Exception as exc:  # Ensembl/UCSC return transient 5xx
            last = exc
            time.sleep(5 + 5 * i)
    raise RuntimeError(f"{url}: {last}")


def gene_models(genes=PANEL) -> dict:
    out = {}
    for g in genes:
        s = _get("https://api.genome.ucsc.edu/search?" + urllib.parse.urlencode({"search": g, "genome": "hg38"}))
        pos = next(x["position"] for m in s["positionMatches"] if m.get("trackName") in ("knownGene", "ncbiRefSeq", "mane")
                   for x in m["matches"] if x.get("posName", "").startswith(g) or g in x.get("description", ""))
        c, se = pos.split(":")
        a, b = (int(v) for v in se.split("-"))
        t = _get("https://api.genome.ucsc.edu/getData/track?" + urllib.parse.urlencode(
            {"genome": "hg38", "track": "mane", "chrom": c, "start": a, "end": b}))
        x = next(r for r in t["mane"] if g in json.dumps(r))
        st = [int(v) for v in x["chromStarts"].strip(",").split(",")]
        sz = [int(v) for v in x["blockSizes"].strip(",").split(",")]
        cs = x["chromStart"]
        out[g] = {"chrom": c.replace("chr", ""), "start": cs + 1, "end": x["chromEnd"], "strand": x["strand"],
                  "tx": x["name"], "exons": [(cs + p + 1, cs + p + q) for p, q in zip(st, sz)]}
    return out


def where(v: dict, pos: int) -> str | None:
    for i, (a, b) in enumerate(v["exons"]):
        if a - 20 <= pos <= b + 20:
            return f"exon{i + 1}" + ("" if a <= pos <= b else "_flank")
    up = v["start"] - 2000 <= pos < v["start"] if v["strand"] == "+" else v["end"] < pos <= v["end"] + 2000
    if up:
        return "promoter2kb"
    return "intron" if v["start"] <= pos <= v["end"] else None


def scan(vcf: str, models: dict):
    byc: dict[str, list] = {}
    for g, v in models.items():
        byc.setdefault(v["chrom"], []).append((g, v))
    keep, counts = [], {g: {"intron": 0, "exonic_or_flank_or_promoter": 0} for g in models}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            c = line.split("\t", 1)[0]
            if c not in byc:
                continue
            p = line.rstrip("\n").split("\t")
            pos = int(p[1])
            for g, v in byc[c]:
                w = where(v, pos)
                if not w:
                    continue
                if w == "intron":
                    counts[g]["intron"] += 1
                    continue
                counts[g]["exonic_or_flank_or_promoter"] += 1
                fmt = dict(zip(p[8].split(":"), p[9].split(":")))
                keep.append({"gene": g, "chrom": c, "pos": pos, "ref": p[3], "alt": p[4], "filter": p[6], "where": w,
                             **{k: fmt.get(k) for k in ("GT", "AD", "DP", "GQ", "PGT", "PID")}})
    return keep, counts


def annotate(keep: list) -> list:
    ann = []
    for i in range(0, len(keep), 50):
        lines = [f"{r['chrom']} {r['pos']} . {r['ref']} {alt} . . ." for r in keep[i:i + 50] for alt in r["alt"].split(",")]
        ann += _get(VEP, json.dumps({"variants": lines}).encode())
    return ann


def summarise(keep: list, ann: list) -> list:
    idx = {(r["chrom"], r["pos"]): r for r in keep}
    out = []
    for a in ann:
        t = next((x for x in a.get("transcript_consequences", []) if x.get("canonical")), {})
        afs = [d[k] for cv in a.get("colocated_variants", []) or [] for d in (cv.get("frequencies") or {}).values()
               for k in ("gnomadg", "gnomade") if k in d]
        af = max(afs) if afs else None
        r = idx.get((a["seq_region_name"], a["start"])) or idx.get((a["seq_region_name"], a["start"] - 1)) or {}
        sai = t.get("spliceai") or {}
        ds = max([v for k, v in sai.items() if k.startswith("DS_")] or [0])
        rare = af is None or af < 0.01
        relevant = t.get("impact") in ("HIGH", "MODERATE") or "splice" in ",".join(t.get("consequence_terms", [])) or ds > 0.2
        if rare and relevant:
            am = t.get("alphamissense") or {}
            out.append({"gene": t.get("gene_symbol"), "hgvsc": t.get("hgvsc"), "hgvsp": t.get("hgvsp"),
                        "consequence": ",".join(t.get("consequence_terms", [])), "gnomad_max_af": af,
                        "alphamissense": am.get("am_pathogenicity"), "revel": t.get("revel"), "cadd": t.get("cadd_phred"),
                        "spliceai_max": ds, "GT": r.get("GT"), "filter": r.get("filter"), "phase_set": r.get("PID"),
                        "clinvar": [cv.get("clin_sig") for cv in a.get("colocated_variants", []) if cv.get("clin_sig")]})
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--workdir", required=True, help="secure directory outside the repository")
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    models = gene_models()
    keep, counts = scan(a.vcf, models)
    ann = annotate(keep)
    json.dump({"records": keep, "vep": ann}, open(os.path.join(a.workdir, "panel_records.json"), "w"))
    rows = summarise(keep, ann)
    print(json.dumps({"panel_genes": len(models), "records_annotated": len(keep), "counts": counts}, indent=1))
    for r in rows:
        print(json.dumps(r))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
