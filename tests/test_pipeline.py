import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from clinvar_landscape import classify  # noqa: E402
from common import load_yaml  # noqa: E402
from score import fit, rank_table  # noqa: E402
from verify_refs import cited_keys  # noqa: E402

WIN = (720, 1050)


def test_classify_truncating_and_missense():
    assert classify("NM_001211.6(BUB1B):c.2364dup (p.Ala789fs)", [], WIN)[0] == "TRUNC"
    assert classify("NM_001211.6(BUB1B):c.100C>T (p.Gln34Ter)", [], WIN)[0] == "TRUNC"
    assert classify("NM_001211.6(BUB1B):c.1058+2T>G", [], WIN)[0] == "TRUNC"
    assert classify("(p.Ile909Thr)", [], WIN) == ("MS_UNSTABLE", 909)
    assert classify("(p.Arg550Gln)", [], WIN) == ("MS_FUNCTIONAL", 550)


def test_every_citation_is_defined():
    refs = load_yaml("references.yaml")
    assert cited_keys(load_yaml("evidence.yaml")) <= set(refs)


def test_dosage_strategy_needs_a_protein_making_allele():
    ev = load_yaml("evidence.yaml")
    s1 = ev["strategies"]["S1_DOSAGE"]
    assert fit(s1, ["TRUNC", "TRUNC"]) == 0
    assert fit(s1, ["TRUNC", "MS_UNSTABLE"]) == 1


def test_primary_ranking_is_deterministic():
    ev, cfg = load_yaml("evidence.yaml"), load_yaml("weights.yaml")
    a, b = rank_table(ev, cfg), rank_table(ev, cfg)
    assert a.equals(b)
    top = a[a.genotype == cfg["primary_genotype"]].sort_values("rank").drug.iloc[0]
    assert top == "niacin"
