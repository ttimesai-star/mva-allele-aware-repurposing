"""Figures for the patient-specific analysis (Submission 2).

figures/patient_alleles.png      BUBR1 map: the proband's two alleles, published MVA1 missense alleles, domains
figures/variant_effect.png       sequence/structure models: N1002K against ClinVar-benign and published unstable alleles
figures/patient_class_sweep.png  rank of each candidate vs P(missense allele is unstable)
"""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import FIGURES, OUTPUTS, PROCESSED, load_yaml  # noqa: E402

INK, MUTED, RED, BLUE, ORANGE, GREEN = "#1f2933", "#6b7785", "#c0392b", "#2471a3", "#e67e22", "#1b7f5a"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})


def patient_alleles():
    pg = load_yaml("patient_genotype.yaml")
    dom = pd.read_csv(PROCESSED / "bubr1_domains.csv")
    lit = pd.read_csv(PROCESSED.parent / "literature_mva1_missense.csv")
    fig, ax = plt.subplots(figsize=(12, 3.6))
    L = 1050
    ax.add_patch(plt.Rectangle((1, -0.12), L - 1, 0.24, fc="#e6edf4", ec=MUTED))
    shown = {"BUB1 N-terminal": ("TPR (BUB1 N-terminal)", "#9fb3c8"), "Protein kinase": ("Pseudokinase domain", "#7fb3d5")}
    for _, r in dom[dom.type == "Domain"].iterrows():
        lab, c = shown.get(r.description, (r.description, "#bcccdc"))
        ax.add_patch(plt.Rectangle((r.start, -0.12), r.end - r.start, 0.24, fc=c, ec=MUTED))
        ax.text((r.start + r.end) / 2, 0, lab, ha="center", va="center", fontsize=9, color=INK)
    for _, r in dom[dom.type == "Motif"].iterrows():
        ax.add_patch(plt.Rectangle((r.start, -0.12), max(6, r.end - r.start), 0.24, fc="#f0b429", ec=MUTED))
        ax.text(r.start, -0.3, r.description, ha="center", va="top", fontsize=7.5, color=MUTED)
    for _, r in lit.iterrows():
        pos = int("".join(ch for ch in r.variant if ch.isdigit()))
        col = BLUE if r.observed_abundance == "low" else "#9fb3c8"
        ax.plot([pos, pos], [-0.12, -0.55], color=col, lw=1)
        ax.plot(pos, -0.58, "o", color=col, ms=6)
        ax.text(pos, -0.68, r.variant, ha="center", va="top", fontsize=7.5, rotation=90, color=col)
    for a in pg["alleles"]:
        pos, col = a["protein_position"], (RED if a["class"].startswith("TRUNC") else BLUE)
        ax.plot([pos, pos], [0.12, 0.62], color=col, lw=2.2)
        ax.plot(pos, 0.66, "v" if a["class"].startswith("TRUNC") else "D", color=col, ms=11)
        ha = "right" if a["class"].startswith("TRUNC") else "left"
        ax.text(pos + (-10 if ha == "right" else 10), 0.8, a['hgvs_p'] + "\n" + a['class_label'], ha=ha, va="bottom", fontsize=9, color=col, weight="bold")
    ax.text(5, 1.2, "Proband (above the bar) vs published MVA1 missense alleles (below; dark = low abundance in cells)",
            fontsize=9.5, color=INK)
    ax.set_xlim(-10, L + 160)
    ax.set_ylim(-1.25, 1.35)
    ax.set_yticks([])
    ax.set_xlabel("BUBR1 residue (UniProt O60566)")
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    fig.savefig(FIGURES / "patient_alleles.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def variant_effect():
    df = pd.read_csv(OUTPUTS / "variant_effect_scores.csv")
    cols = [c for c in ("esm1v_llr", "esm2_llr", "thermompnn_ddg") if c in df and df[c].notna().any()]
    titles = {"esm1v_llr": "ESM-1v log-likelihood ratio\n(lower = less tolerated)",
              "esm2_llr": "ESM-2 log-likelihood ratio\n(lower = less tolerated)",
              "thermompnn_ddg": "ThermoMPNN predicted ddG, kcal/mol\n(higher = more destabilising)"}
    groups = [("clinvar_benign", "ClinVar benign /\nlikely benign", "#9fb3c8"),
              ("lit_low", "Published MVA1,\nlow abundance", BLUE), ("lit_normal", "Published MVA1,\nnormal abundance", "#7fb3d5")]
    fig, axes = plt.subplots(1, len(cols), figsize=(4.3 * len(cols), 4.2))
    axes = np.atleast_1d(axes)
    rng = np.random.default_rng(1)
    for ax, c in zip(axes, cols):
        for k, (g, lab, col) in enumerate(groups):
            v = df[(df.set == g)][c].dropna()
            ax.scatter(k + rng.uniform(-0.12, 0.12, len(v)), v, color=col, s=26, alpha=0.85, edgecolor="white", lw=0.5)
            if g in ("lit_low", "lit_normal"):
                sub = df[(df.set == g)].dropna(subset=[c]).sort_values(c)
                span = ax.get_ylim() if False else (df[c].min(), df[c].max())
                gap = (span[1] - span[0]) * 0.055
                last = -1e9
                for _, r in sub.iterrows():
                    yv = max(r[c], last + gap)
                    ax.annotate(r.variant, (k, r[c]), (k + 0.2, yv), fontsize=7, va="center", color=MUTED,
                                arrowprops=dict(arrowstyle="-", color="#cbd2d9", lw=0.6))
                    last = yv
        q = df[df.set == "query"][c].dropna()
        ax.axhline(q.iloc[0], color=RED, lw=1.6, ls="--")
        ax.text(-0.45, q.iloc[0], " N1002K (proband)", color=RED, fontsize=8.5, va="bottom", ha="left", weight="bold")
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([g[1] for g in groups], fontsize=8)
        ax.set_title(titles[c], fontsize=9.5)
        for s in ("right", "top"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURES / "variant_effect.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def class_sweep():
    sw = pd.read_csv(OUTPUTS / "patient_class_sweep.csv")
    p0 = load_yaml("weights.yaml")["patient"]["p_missense_unstable"]
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    top = sw[sw.p_missense_unstable == p0].sort_values("rank").drug.tolist()[:6]
    pal = [GREEN, BLUE, "#8e44ad", ORANGE, "#2e9e74", "#b9770e"]
    for d, c in zip(top, pal):
        s = sw[sw.drug == d]
        ax.plot(s.p_missense_unstable, s["score"], color=c, lw=2, label=d)
    ax.axvline(p0, color=MUTED, ls=":", lw=1.2)
    ax.text(p0, ax.get_ylim()[1], f" prior used = {p0}", color=MUTED, fontsize=8.5, va="top")
    ax.set_xlabel("P(the proband's missense allele makes unstable, fast-turnover protein)")
    ax.set_ylabel("expected allele-aware score")
    ax.legend(frameon=False, fontsize=8.5, ncol=2)
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)
    fig.savefig(FIGURES / "patient_class_sweep.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    patient_alleles()
    if (OUTPUTS / "variant_effect_scores.csv").exists():
        variant_effect()
    class_sweep()
    return 0


if __name__ == "__main__":
    sys.exit(main())
