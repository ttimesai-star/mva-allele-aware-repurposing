"""Allele-aware ranking of candidate medicines + weight/score sensitivity analysis.

score(drug, genotype) = base(drug) * (fit_floor + (1 - fit_floor) * fit(strategy, genotype))
base = sum_i w_i * s_i / 3  -  w_risk * risk / 3          (i over mech, model, human, ped, test)
fit(strategy, genotype) = max over the two alleles of strategy.allele_fit[allele]
        (one responsive allele is enough for a dosage-type strategy)

Outputs
  outputs/ranking_by_genotype.csv
  outputs/rank_stability_<primary>.csv   (P(rank 1), P(top 3), median rank over n_draws)
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from common import OUTPUTS, load_yaml

DIMS = ["mech", "model", "human", "ped", "test"]


def fit(strategy: dict, genotype: list[str]) -> float:
    return max(float(strategy["allele_fit"].get(a, 0.0)) for a in genotype)


def base_score(scores: dict, w: dict) -> float:
    return sum(w[d] * scores[d] / 3 for d in DIMS) - w["risk"] * scores["risk"] / 3


def rank_table(ev: dict, cfg: dict) -> pd.DataFrame:
    rows = []
    for gname, gt in cfg["genotypes"].items():
        for c in ev["candidates"]:
            st = ev["strategies"][c["strategy"]]
            f = fit(st, gt)
            b = base_score(c["scores"], cfg["weights"])
            rows.append({"genotype": gname, "drug": c["drug"], "strategy": c["strategy"], "base": round(b, 3),
                         "fit": f, "score": round(b * (cfg["fit_floor"] + (1 - cfg["fit_floor"]) * f), 3)})
    df = pd.DataFrame(rows)
    df["rank"] = df.groupby("genotype")["score"].rank(ascending=False, method="min").astype(int)
    return df.sort_values(["genotype", "rank"])


def stability(ev: dict, cfg: dict, genotype: str) -> pd.DataFrame:
    rng = np.random.default_rng(cfg["seed"])
    gt = cfg["genotypes"][genotype]
    keys = DIMS + ["risk"]
    w0 = np.array([cfg["weights"][k] for k in keys])
    cands = ev["candidates"]
    S = np.array([[c["scores"][k] for k in keys] for c in cands], dtype=float)
    F = np.array([fit(ev["strategies"][c["strategy"]], gt) for c in cands])
    ranks = np.zeros((cfg["n_draws"], len(cands)), dtype=int)
    for t in range(cfg["n_draws"]):
        w = rng.dirichlet(w0 / w0.sum() * cfg["dirichlet_concentration"]) * w0.sum()
        s = np.clip(S + rng.uniform(-cfg["score_jitter"], cfg["score_jitter"], S.shape), 0, 3)
        base = (s[:, :5] / 3) @ w[:5] - w[5] * s[:, 5] / 3
        sc = base * (cfg["fit_floor"] + (1 - cfg["fit_floor"]) * F)
        order = np.argsort(-sc)
        ranks[t, order] = np.arange(1, len(cands) + 1)
    return pd.DataFrame({"drug": [c["drug"] for c in cands],
                         "p_rank1": (ranks == 1).mean(0).round(3),
                         "p_top3": (ranks <= 3).mean(0).round(3),
                         "median_rank": np.median(ranks, 0)}).sort_values("p_top3", ascending=False)


# ---------------------------------------------------------------------------------------------
# Patient-specific layer (Submission 2): the confirmed genotype of the hackathon proband, uncertainty
# about the functional class of the missense allele, and phenotype-driven risk modifiers.
# ---------------------------------------------------------------------------------------------
def _risk_mod(ev: dict) -> dict:
    out: dict[str, float] = {}
    for m in ev.get("patient_context_modifiers", []):
        out[m["drug"]] = out.get(m["drug"], 0.0) + float(m["delta_risk"])
    return out


def patient_scores(ev: dict, cfg: dict, p_unstable: float | None = None, use_modifiers: bool = True) -> pd.DataFrame:
    pt = cfg["patient"]
    p = pt["p_missense_unstable"] if p_unstable is None else p_unstable
    mod = _risk_mod(ev) if use_modifiers else {}
    rows = []
    for c in ev["candidates"]:
        sc = dict(c["scores"])
        sc["risk"] = min(3.0, sc["risk"] + mod.get(c["drug"], 0.0))
        b = base_score(sc, cfg["weights"])
        st = ev["strategies"][c["strategy"]]
        f1, f0 = fit(st, pt["alleles"]), fit(st, pt["alleles_if_missense_stable"])
        s1 = b * (cfg["fit_floor"] + (1 - cfg["fit_floor"]) * f1)
        s0 = b * (cfg["fit_floor"] + (1 - cfg["fit_floor"]) * f0)
        rows.append({"drug": c["drug"], "strategy": c["strategy"], "risk_used": sc["risk"],
                     "score_if_unstable": round(s1, 3), "score_if_stable": round(s0, 3),
                     "score": round(p * s1 + (1 - p) * s0, 3)})
    df = pd.DataFrame(rows)
    df["rank"] = df["score"].rank(ascending=False, method="min").astype(int)
    return df.sort_values("rank")


def patient_stability(ev: dict, cfg: dict) -> pd.DataFrame:
    """Monte Carlo: perturbed weights and scores (as `stability`) plus, in every draw, the missense allele is
    drawn as unstable with probability p_missense_unstable and as stable otherwise."""
    rng = np.random.default_rng(cfg["seed"] + 1)
    pt = cfg["patient"]
    keys = DIMS + ["risk"]
    mod = _risk_mod(ev)
    w0 = np.array([cfg["weights"][k] for k in keys])
    cands = ev["candidates"]
    S = np.array([[c["scores"][k] + (mod.get(c["drug"], 0.0) if k == "risk" else 0) for k in keys] for c in cands], float)
    F1 = np.array([fit(ev["strategies"][c["strategy"]], pt["alleles"]) for c in cands])
    F0 = np.array([fit(ev["strategies"][c["strategy"]], pt["alleles_if_missense_stable"]) for c in cands])
    ranks = np.zeros((cfg["n_draws"], len(cands)), dtype=int)
    for t in range(cfg["n_draws"]):
        w = rng.dirichlet(w0 / w0.sum() * cfg["dirichlet_concentration"]) * w0.sum()
        s = np.clip(S + rng.uniform(-cfg["score_jitter"], cfg["score_jitter"], S.shape), 0, 3)
        base = (s[:, :5] / 3) @ w[:5] - w[5] * s[:, 5] / 3
        F = F1 if rng.random() < pt["p_missense_unstable"] else F0
        sc = base * (cfg["fit_floor"] + (1 - cfg["fit_floor"]) * F)
        order = np.argsort(-sc)
        ranks[t, order] = np.arange(1, len(cands) + 1)
    return pd.DataFrame({"drug": [c["drug"] for c in cands],
                         "p_rank1": (ranks == 1).mean(0).round(3),
                         "p_top3": (ranks <= 3).mean(0).round(3),
                         "median_rank": np.median(ranks, 0)}).sort_values("p_top3", ascending=False)


def class_sweep(ev: dict, cfg: dict) -> pd.DataFrame:
    """Rank of every candidate as a function of P(missense allele is unstable)."""
    rows = []
    for p in np.round(np.linspace(0, 1, 21), 2):
        d = patient_scores(ev, cfg, p_unstable=float(p))
        for _, r in d.iterrows():
            rows.append({"p_missense_unstable": p, "drug": r["drug"], "score": r["score"], "rank": r["rank"]})
    return pd.DataFrame(rows)


def main() -> int:
    ev, cfg = load_yaml("evidence.yaml"), load_yaml("weights.yaml")
    df = rank_table(ev, cfg)
    df.to_csv(OUTPUTS / "ranking_by_genotype.csv", index=False)
    g = cfg["primary_genotype"]
    print(df[df.genotype == g].to_string(index=False))
    st = stability(ev, cfg, g)
    st.to_csv(OUTPUTS / f"rank_stability_{g}.csv", index=False)
    print(st.to_string(index=False))
    if "patient" in cfg:
        ps = patient_scores(ev, cfg)
        ps.to_csv(OUTPUTS / "patient_ranking.csv", index=False)
        patient_scores(ev, cfg, use_modifiers=False).to_csv(OUTPUTS / "patient_ranking_no_context.csv", index=False)
        pst = patient_stability(ev, cfg)
        pst.to_csv(OUTPUTS / "patient_rank_stability.csv", index=False)
        sw = class_sweep(ev, cfg)
        sw.to_csv(OUTPUTS / "patient_class_sweep.csv", index=False)
        print("\nPATIENT ranking (genotype + context modifiers):")
        print(ps.to_string(index=False))
        print(pst.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
