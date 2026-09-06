import numpy as np
import pandas as pd
import pytest

from contextshift import stats
from contextshift.stats import GLOBAL, PER_PAIR, add_qvalues, benjamini_hochberg


def test_bh_matches_known_values():
    p = np.array([0.01, 0.02, 0.03, 0.04, 0.05])
    q = benjamini_hochberg(p)
    assert q == pytest.approx([0.05, 0.05, 0.05, 0.05, 0.05])


def test_bh_is_monotone():
    rng = np.random.default_rng(0)
    q = benjamini_hochberg(rng.uniform(size=500))
    assert np.all(np.diff(np.sort(q)) >= -1e-12)


def test_bh_never_below_pvalue():
    p = np.array([0.001, 0.2, 0.7])
    assert np.all(benjamini_hochberg(p) >= p - 1e-12)


def test_bh_handles_nan():
    q = benjamini_hochberg(np.array([0.01, np.nan, 0.02]))
    assert np.isnan(q[1])
    assert np.isfinite(q[0]) and np.isfinite(q[2])


def sites_frame(n_pairs=2, n_sites=50):
    rows = []
    for pair in range(n_pairs):
        for c in range(n_sites):
            rows.append(
                {
                    "family": "F",
                    "partition": "context",
                    "group_a": f"A{pair}",
                    "group_b": f"B{pair}",
                    "column": c,
                    "test": "type2",
                    "posterior": None,
                    "pvalue": (c + 1) / (n_sites * 20),
                }
            )
    return pd.DataFrame(rows)


def test_global_correction_is_stricter_than_per_pair():
    df = sites_frame()
    g = add_qvalues(df, scope=GLOBAL)["qvalue"].to_numpy()
    p = add_qvalues(df, scope=PER_PAIR)["qvalue"].to_numpy()
    assert np.nanmean(g) >= np.nanmean(p)


def test_unknown_scope_rejected():
    with pytest.raises(ValueError):
        add_qvalues(sites_frame(), scope="whatever")


def test_testing_burden_counts_every_comparison():
    burden = stats.testing_burden(sites_frame(n_pairs=3, n_sites=10))
    assert len(burden) == 3
    assert set(burden["n_tests"]) == {10}
