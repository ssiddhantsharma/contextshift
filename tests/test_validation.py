"""Properties the pipeline must have, checked against method definitions."""

from pathlib import Path

import pytest

from contextshift.stages import diverge

DATA = Path(__file__).parent / "data" / "diverge_casp"
ALN, T1, T2, T3 = DATA / "CASP.aln", DATA / "cl1.tree", DATA / "cl2.tree", DATA / "cl3.tree"

pytestmark = pytest.mark.skipif(not diverge.available(), reason="diverge not installed")


def gu99(*trees, names):
    d = diverge._import_diverge()
    diverge.apply_shim()
    return d.Gu99(str(ALN), *[str(t) for t in trees], cluster_name=names)


def test_pairwise_theta_matches_a_multi_cluster_run():
    """Running one pair alone must equal extracting it from a larger run.

    This library loops over group pairs rather than handing DIVERGE every
    cluster at once; that is only sound if the pairwise estimate is unchanged.
    """
    two = gu99(T1, T2, names=["A", "B"])
    three = gu99(T1, T2, T3, names=["A", "B", "C"])
    assert two.summary.loc["ThetaML", "A/B"] == three.summary.loc["ThetaML", "A/B"]


def test_pairwise_per_site_posteriors_match_a_multi_cluster_run():
    two = gu99(T1, T2, names=["A", "B"])
    three = gu99(T1, T2, T3, names=["A", "B", "C"])
    a, b = two.results["A/B"], three.results["A/B"]
    assert list(a.index) == list(b.index)
    assert float((a - b).abs().max()) == 0.0


def test_kept_positions_do_not_depend_on_the_other_clusters():
    two = gu99(T1, T2, names=["A", "B"])
    three = gu99(T1, T2, T3, names=["A", "B", "C"])
    assert list(two.results.index) == list(three.results.index)
