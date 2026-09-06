"""Runs DIVERGE for real, against the fixtures shipped with DIVERGE itself.

Every constant asserted here was measured from a live 4.1.0 run, not taken
from documentation. If a future DIVERGE changes any of them, this fails loudly
rather than letting a silently different convention through.
"""

from pathlib import Path

import pytest

from contextshift.partition import Partition, Provenance
from contextshift.schema import COMPARISONS, SITES
from contextshift.stages import diverge

pytestmark = pytest.mark.skipif(not diverge.available(), reason="diverge not installed")

DATA = Path(__file__).parent / "data" / "diverge_casp"
ALN, T1, T2 = DATA / "CASP.aln", DATA / "cl1.tree", DATA / "cl2.tree"

N_KEPT = 781
FIRST_POSITION = 132
LAST_POSITION = 1460


def run():
    return diverge.run_pair(ALN, T1, T2, "CASP", "context", "A", "B")


def test_runs_and_conforms_to_schema():
    sites, comparisons, _ = run()
    SITES.validate(sites)
    COMPARISONS.validate(comparisons)


def test_positions_are_sparse_not_a_dense_range():
    sites, _, _ = run()
    t2 = sites[sites["test"] == diverge.TYPE2]
    assert len(t2) == N_KEPT
    assert int(t2["column"].min()) == FIRST_POSITION
    assert int(t2["column"].max()) == LAST_POSITION
    # the whole point: DIVERGE drops columns, so this is not 0..N-1
    assert len(t2) < LAST_POSITION - FIRST_POSITION + 1


def test_value_is_a_posterior_not_a_pvalue():
    sites, _, _ = run()
    post = sites["posterior"].dropna()
    assert post.min() >= 0.0
    assert post.max() <= 1.0
    # a p-value column must not be invented where the method produces none
    assert "pvalue" not in sites.columns


def test_type1_and_type2_both_returned_on_the_same_positions():
    sites, _, _ = run()
    t1 = sites[sites["test"] == diverge.TYPE1]
    t2 = sites[sites["test"] == diverge.TYPE2]
    assert len(t1) == len(t2) == N_KEPT
    assert list(t1["column"]) == list(t2["column"])


def test_shim_use_is_recorded_when_needed():
    _, _, notes = run()
    if diverge.needs_shim():
        assert any("shim" in n for n in notes)


def test_comparison_coefficients_are_per_pair_not_per_site():
    _, comparisons, _ = run()
    params = set(comparisons["parameter"])
    assert "AlphaML" in params
    assert {"Theta-II", "ThetaML"} & params
    # one row per parameter per test, never one per alignment column
    assert len(comparisons) < 100


def test_group_labels_survive_into_the_column_name():
    sites, _, _ = run()
    assert set(sites["group_a"]) == {"A"}
    assert set(sites["group_b"]) == {"B"}


def test_run_partition_reports_skips_and_notes():
    p = Partition(
        name="context",
        labels={"A": "A", "B": "B", "C": "C"},
        provenance=Provenance(source="test"),
    )
    sites, comparisons, skipped, notes = diverge.run_partition(
        ALN, {"A": T1, "B": T2}, p, family="CASP"
    )
    assert len(sites) == 2 * N_KEPT
    # the A|C and B|C pairs have no tree and must be reported, not dropped
    assert len(skipped) == 2
    assert all("missing conformant tree" in s for s in skipped)
