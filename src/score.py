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


def main() -> int:
    ev, cfg = load_yaml("evidence.yaml"), load_yaml("weights.yaml")
    df = rank_table(ev, cfg)
    df.to_csv(OUTPUTS / "ranking_by_genotype.csv", index=False)
    g = cfg["primary_genotype"]
    print(df[df.genotype == g].to_string(index=False))
    st = stability(ev, cfg, g)
    st.to_csv(OUTPUTS / f"rank_stability_{g}.csv", index=False)
    print(st.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
