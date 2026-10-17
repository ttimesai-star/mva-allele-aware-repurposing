"""Figures for the report and the pitch video (PNG, 200 dpi)."""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

from common import FIGURES, OUTPUTS, PROCESSED, load_yaml  # noqa: E402

INK = "#1f2933"
MUTED = "#6b7785"
C = {"TRUNC": "#c0392b", "TRUNC_NONSENSE": "#e67e22", "MS_UNSTABLE": "#2471a3", "MS_FUNCTIONAL": "#7fb3d5",
     "REG": "#8e44ad", "OTHER": "#95a5a6"}
STRAT_C = {"S1_DOSAGE": "#1b7f5a", "S1b_PROTEOSTASIS": "#2e9e74", "S1c_READTHROUGH": "#7cc6a4",
           "S2_MTOR": "#2471a3", "S2b_IL1": "#5dade2", "S3_SENOLYTIC": "#8e44ad", "S4_ANEU_SELECTIVE": "#b9770e"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})


def box(ax, x, y, w, h, text, fc, ec=None, fs=9.5, color="white", weight="bold"):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.012,rounding_size=0.018",
                                fc=fc, ec=ec or fc, lw=1.2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=color, weight=weight, wrap=True)


def arrow(ax, p, q, color=MUTED, style="-|>", lw=1.4, ls="-"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=12, color=color, lw=lw, linestyle=ls,
                                 shrinkA=4, shrinkB=4))


def mechanism_map():
    fig, ax = plt.subplots(figsize=(13, 7.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.01, 1.06)
    ax.axis("off")
    dis = "#34495e"
    # disease chain (left -> right)
    chain = [(0.09, 0.62, "Biallelic BUB1B\ntruncating + missense"),
             (0.27, 0.62, "Low BUBR1\n(no transcript +\nunstable missense)"),
             (0.45, 0.62, "Weak spindle\ncheckpoint +\nattachment errors"),
             (0.63, 0.62, "Chromosome\nmissegregation\n(mosaic aneuploidy)")]
    for x, y, t in chain:
        box(ax, x, y, 0.15, 0.15, t, dis)
    for (x1, y1, _), (x2, y2, _) in zip(chain, chain[1:]):
        arrow(ax, (x1 + 0.075, y1), (x2 - 0.075, y2), color=dis, lw=2)
    outs = [(0.86, 0.86, "Embryonal cancers\n(RMS, Wilms, leukaemia)", "#922b21"),
            (0.86, 0.62, "Senescence + SASP\n(p53-p21, p16)", "#6c3483"),
            (0.86, 0.38, "Progeroid decline\n(lens, fat, muscle)", "#6c3483")]
    for x, y, t, c in outs:
        box(ax, x, y, 0.18, 0.12, t, c)
    arrow(ax, (0.705, 0.66), (0.77, 0.84), color=dis, lw=2)
    arrow(ax, (0.705, 0.62), (0.77, 0.62), color=dis, lw=2)
    arrow(ax, (0.86, 0.56), (0.86, 0.44), color="#6c3483", lw=2)
    box(ax, 0.60, 0.40, 0.2, 0.08, "mTORC1 hyperactivity in muscle\n(truncating allele, mice)", "#5d6d7e", fs=8.3)
    arrow(ax, (0.12, 0.545), (0.50, 0.41), color="#5d6d7e", ls="--")
    arrow(ax, (0.70, 0.40), (0.77, 0.385), color="#5d6d7e", ls="--")
    iv = [(0.155, 0.20, "S1  NAD+ -> SIRT2\nstabilise BUBR1\nniacin", STRAT_C["S1_DOSAGE"], (0.25, 0.545)),
          (0.385, 0.20, "S1b  heat-shock response\nhelp missense fold\narimoclomol", STRAT_C["S1b_PROTEOSTASIS"], (0.29, 0.545)),
          (0.60, 0.20, "S2  dampen mTORC1 / SASP\nlow-dose intermittent\neverolimus", STRAT_C["S2_MTOR"], (0.60, 0.36)),
          (0.84, 0.17, "S3  clear senescent cells\nhit-and-run\ndasatinib (+ quercetin)", STRAT_C["S3_SENOLYTIC"], (0.955, 0.60))]
    for x, y, t, c, tgt in iv:
        box(ax, x, y, 0.18, 0.13, t, c, fs=8.4)
        src = (x + 0.085, y + 0.03) if x > 0.8 else (x, y + 0.065)
        arrow(ax, src, tgt, color=c, style="-[", lw=1.8)
    box(ax, 0.30, 0.04, 0.56, 0.06, "Do-no-harm layer: tubulin agents, HSP90 inhibitors, radiation  ->  flag and monitor",
        "#fdebd0", ec="#b9770e", color="#7e5109", fs=8.6)
    ax.text(0.01, 1.04, "Mechanism map: where each repurposing hypothesis acts", fontsize=13, weight="bold", color=INK, va="top")
    ax.text(0.01, 1.0, "solid arrows = shown in BubR1-insufficient mice or MVA patient cells; dashed = shown in one mouse genotype; bars = proposed drug action",
            fontsize=9, color=MUTED, va="top")
    fig.savefig(FIGURES / "mechanism_map.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def allele_lollipop():
    dom = pd.read_csv(PROCESSED / "bubr1_domains.csv")
    v = pd.read_csv(PROCESSED / "clinvar_mva_variants.csv")
    v = v[(v.gene == "BUB1B") & v.protein_pos.notna()]
    lit = pd.read_csv(OUTPUTS / "kinase_window_rule_check.csv")
    ev = load_yaml("evidence.yaml")
    fig, ax = plt.subplots(figsize=(12, 4.2))
    L = 1050
    ax.add_patch(plt.Rectangle((1, -0.12), L, 0.24, fc="#d5dbdb", ec="none"))
    for _, d in dom[dom.type == "Domain"].iterrows():
        ax.add_patch(plt.Rectangle((d.start, -0.18), d.end - d.start, 0.36, fc="#5d6d7e", ec="none"))
        ax.text((d.start + d.end) / 2, 0, d.description, ha="center", va="center", color="white", fontsize=8.5, weight="bold")
    w = ev["kinase_window"]
    ax.axvspan(w["start"], w["end"], ymin=0.0, ymax=1, color="#2471a3", alpha=0.07)
    ax.text(w["start"] + 5, -1.12, "'in or near kinase domain' window -> unstable missense", color="#2471a3", fontsize=8.5)
    for i, (_, r) in enumerate(v.sort_values("protein_pos").iterrows()):
        h = 0.45 + 0.55 * ((i % 4) / 3)
        ax.plot([r.protein_pos] * 2, [0.18, h], color=C[r.allele_class], lw=0.8)
        ax.scatter(r.protein_pos, h, s=22, color=C[r.allele_class], zorder=3)
    for k, (_, r) in enumerate(lit.iterrows()):
        pos = int("".join(ch for ch in r.variant if ch.isdigit()))
        col = C["MS_UNSTABLE"] if r.observed_abundance == "low" else C["MS_FUNCTIONAL"]
        ax.plot([pos] * 2, [-0.18, -0.7], color=col, lw=1)
        ax.scatter(pos, -0.7, s=40, color=col, marker="D", zorder=3, edgecolor="white")
        ax.text(pos, -0.84 - 0.13 * (k % 2), r.variant, ha="center", fontsize=7.5, color=INK, rotation=0)
    ax.set_xlim(0, L + 10)
    ax.set_ylim(-1.2, 1.2)
    ax.set_yticks([])
    ax.set_xlabel("BUBR1 residue (UniProt O60566)")
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    n = v.allele_class.value_counts()
    ax.set_title(f"Above: ClinVar pathogenic/likely pathogenic BUB1B variants with a protein position (n={len(v)}; "
                 f"truncating {n.get('TRUNC', 0)}, nonsense {n.get('TRUNC_NONSENSE', 0)}, missense {n.get('MS_FUNCTIONAL', 0) + n.get('MS_UNSTABLE', 0)})\n"
                 "Below: published MVA1 missense alleles, coloured by measured protein abundance (Suijkerbuijk 2010)",
                 fontsize=10, color=INK, loc="left")
    handles = [plt.Line2D([], [], marker="o", ls="", color=C[k], label=lab) for k, lab in
               [("TRUNC", "frameshift / splice"), ("TRUNC_NONSENSE", "nonsense"), ("MS_UNSTABLE", "missense, low protein"),
                ("MS_FUNCTIONAL", "missense, normal protein")]]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=8.5, ncol=4, bbox_to_anchor=(1, 1.02))
    fig.savefig(FIGURES / "bubr1_alleles.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def ranking_heatmap():
    r = pd.read_csv(OUTPUTS / "ranking_by_genotype.csv")
    cfg = load_yaml("weights.yaml")
    piv = r.pivot(index="drug", columns="genotype", values="score")
    order = piv[cfg["primary_genotype"]].sort_values(ascending=False).index
    piv = piv.loc[order, list(cfg["genotypes"])]
    labels = {"G1_TRUNC_x_MSUNSTABLE": "truncating +\nunstable missense\n(most common)",
              "G2_NONSENSE_x_MSUNSTABLE": "nonsense +\nunstable missense",
              "G3_TRUNC_x_MSFUNCTIONAL": "truncating +\nstable missense",
              "G4_TRUNC_x_REG": "truncating +\nlow-expression", "G5_MSUNSTABLE_x_MSUNSTABLE": "two unstable\nmissense"}
    fig, ax = plt.subplots(figsize=(10, 5.6))
    im = ax.imshow(piv.values, cmap="YlGnBu", vmin=0, vmax=0.7, aspect="auto")
    ax.set_xticks(range(piv.shape[1]), [labels[c] for c in piv.columns], fontsize=8.8)
    strat = r.drop_duplicates("drug").set_index("drug")["strategy"]
    ax.set_yticks(range(piv.shape[0]), [f"{d}  [{strat[d].split('_')[0]}]" for d in piv.index])
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            val = piv.values[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8.5, color="white" if val > 0.45 else INK)
    ax.set_title("Allele-aware score by BUB1B genotype (higher = stronger, safer, more testable)", loc="left", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.savefig(FIGURES / "ranking_heatmap.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def stability_bars():
    cfg = load_yaml("weights.yaml")
    st = pd.read_csv(OUTPUTS / f"rank_stability_{cfg['primary_genotype']}.csv").sort_values("p_top3")
    r = pd.read_csv(OUTPUTS / "ranking_by_genotype.csv").drop_duplicates("drug").set_index("drug")["strategy"]
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    ax.barh(st.drug, st.p_top3, color=[STRAT_C[r[d]] for d in st.drug])
    for y, (p, p1) in enumerate(zip(st.p_top3, st.p_rank1)):
        ax.text(p + 0.01, y, f"top-3 {p:.0%}  |  #1 {p1:.0%}", va="center", fontsize=8.5, color=INK)
    ax.set_xlim(0, 1.25)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0], ["0%", "25%", "50%", "75%", "100%"])
    ax.set_xlabel(f"share of {cfg['n_draws']} perturbed rankings (random weights + score noise)")
    ax.set_title("How robust is the ranking? (truncating + unstable missense genotype)", loc="left", fontsize=11)
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)
    fig.savefig(FIGURES / "rank_stability.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    mechanism_map()
    allele_lollipop()
    ranking_heatmap()
    stability_bars()
    print("figures:", ", ".join(p.name for p in sorted(FIGURES.glob("*.png"))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
