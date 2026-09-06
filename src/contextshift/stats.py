"""Multiple-testing correction.

Sites are tested per family, per group pair, per column, so the default
correction is global over the whole table.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .schema import SITES

GLOBAL = "global"
PER_FAMILY = "family"
PER_PAIR = "pair"

_SCOPES = {
    GLOBAL: [],
    PER_FAMILY: ["family"],
    PER_PAIR: ["family", "partition", "group_a", "group_b", "test"],
}


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return out

    vals = p[finite]
    n = vals.size
    order = np.argsort(vals)
    ranked = vals[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    restored = np.empty(n)
    restored[order] = q
    out[finite] = restored
    return out


def add_qvalues(sites: pd.DataFrame, scope: str = GLOBAL) -> pd.DataFrame:
    if scope not in _SCOPES:
        raise ValueError(f"scope must be one of {sorted(_SCOPES)}, got {scope!r}")
    df = SITES.validate(sites.copy())
    keys = _SCOPES[scope]

    if not keys:
        df["qvalue"] = benjamini_hochberg(df["pvalue"].to_numpy())
    else:
        df["qvalue"] = np.nan
        for _, idx in df.groupby(keys, dropna=False).groups.items():
            rows = df.loc[idx]
            df.loc[idx, "qvalue"] = benjamini_hochberg(rows["pvalue"].to_numpy())

    df["qvalue"] = df["qvalue"].astype("float64")
    return SITES.validate(df)


def testing_burden(sites: pd.DataFrame) -> pd.DataFrame:
    """Tests contributed by each comparison."""
    df = SITES.validate(sites.copy())
    return (
        df.groupby(["family", "partition", "group_a", "group_b", "test"], dropna=False)
        .agg(n_tests=("pvalue", "size"), n_testable=("pvalue", "count"))
        .reset_index()
    )
