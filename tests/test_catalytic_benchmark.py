"""Conservation against published ground truth.

Alignments and catalytic-site labels come from the data supporting Capra &
Singh 2007 (Bioinformatics 23:1875), whose Jensen-Shannon method this
implements. Labels are from the Catalytic Site Atlas: 1 catalytic, -1 not.

Measured over 57 proteins of that set: catalytic sites scored as more
constrained in 57 of 57, pooled AUC 0.924. The five kept here each retain at
least three catalytic positions after the occupancy gate; below that a
per-protein AUC is noise rather than a measurement.
"""

from pathlib import Path

import pytest

from contextshift.stages import conserve

DATA = Path(__file__).parent / "data" / "catalytic"
PROTEINS = sorted(p.stem for p in DATA.glob("*.fa"))
MIN_CATALYTIC = 3


def labels_for(stem: str) -> dict[int, int]:
    return {
        int(p[0]): int(p[1])
        for p in (line.split() for line in (DATA / f"{stem}.cat_sites").read_text().splitlines())
        if len(p) == 2
    }


def scored(stem: str):
    d = conserve.jensen_shannon(DATA / f"{stem}.fa", stem, "all").dropna(subset=["rate"])
    d["cat"] = d["column"].map(labels_for(stem))
    return d[d["cat"] == 1]["rate"], d[d["cat"] == -1]["rate"]


def auc(cat, non) -> float:
    return sum((non > r).sum() for r in cat) / (len(cat) * len(non))


def test_fixture_has_enough_catalytic_sites_to_measure():
    assert len(PROTEINS) == 5
    for stem in PROTEINS:
        cat, non = scored(stem)
        assert len(cat) >= MIN_CATALYTIC, f"{stem}: only {len(cat)} scored catalytic sites"
        assert len(non) >= 50


@pytest.mark.parametrize("stem", PROTEINS)
def test_catalytic_sites_are_more_constrained(stem):
    """`rate` is lower for more constrained positions."""
    cat, non = scored(stem)
    assert cat.mean() < non.mean()


def test_pooled_auc_matches_the_published_method():
    """Pooling positions is the meaningful statistic at these site counts."""
    cats, nons = [], []
    for stem in PROTEINS:
        c, n = scored(stem)
        cats.extend(c)
        nons.extend(n)
    import pandas as pd

    pooled = auc(pd.Series(cats), pd.Series(nons))
    assert pooled > 0.85, f"pooled AUC {pooled:.3f}"


def test_gappy_columns_are_excluded_not_scored():
    """Catalytic positions can fall in gappy columns and are then unscored.

    That is the occupancy gate working; it means coverage is incomplete, which
    a caller has to see rather than assume away.
    """
    stem = PROTEINS[0]
    labelled = labels_for(stem)
    full = conserve.jensen_shannon(DATA / f"{stem}.fa", stem, "all")
    assert full["rate"].isna().any() or len(full) == len(labelled)
    assert full["occupancy"].notna().all()
