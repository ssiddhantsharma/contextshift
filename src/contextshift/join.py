"""Join divergence calls to conservation to classify each alignment column.

Divergence alone cannot say whether a site matters; conservation alone cannot
see a between-group shift because it averages over the groups. The cross of the
two is the result:

    conserved everywhere, not divergent   -> core        (catalytic; negative control)
    conserved within groups, differs      -> determinant (the answer)
    constrained in one group only         -> relaxed     (Type-I)
    variable everywhere                   -> variable
    conserved everywhere AND divergent    -> suspect     (usually misalignment)

Columns are never dropped. A column whose inputs are missing is emitted with
status set, so a site cannot vanish by accident.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .schema import CONSERVATION, SITES

CORE = "core"
DETERMINANT = "determinant"
RELAXED = "relaxed"
VARIABLE = "variable"
SUSPECT = "suspect"
UNKNOWN = "unknown"

SCOPE_ALL = "all"


@dataclass(frozen=True)
class JoinThresholds:
    """Rates are relative evolutionary rates: lower means more constrained."""

    conserved_rate: float = 0.5
    variable_rate: float = 1.5
    alpha: float = 0.05
    min_occupancy: float = 0.5


def _significant(df: pd.DataFrame, alpha: float) -> pd.Series:
    col = "qvalue" if "qvalue" in df.columns and df["qvalue"].notna().any() else "pvalue"
    return df[col].notna() & (df[col] <= alpha)


def classify(
    sites: pd.DataFrame,
    conservation: pd.DataFrame,
    thresholds: JoinThresholds | None = None,
) -> pd.DataFrame:
    """Return one row per (family, partition, group pair, column) with a class."""
    t = thresholds or JoinThresholds()
    sites = SITES.validate(sites.copy())
    conservation = CONSERVATION.validate(conservation.copy())

    sites["significant"] = _significant(sites, t.alpha)

    calls = (
        sites.pivot_table(
            index=["family", "partition", "group_a", "group_b", "column"],
            columns="test",
            values="significant",
            aggfunc="max",
        )
        .fillna(False)
        .astype(bool)
        .reset_index()
    )
    for test in ("type1", "type2"):
        if test not in calls.columns:
            calls[test] = False

    cons = conservation.set_index(["family", "column", "scope"])

    def rate_for(family: str, column: int, scope: str) -> float:
        try:
            return float(cons.loc[(family, column, scope), "rate"])
        except (KeyError, TypeError):
            return float("nan")

    def occupancy_for(family: str, column: int, scope: str) -> float:
        if "occupancy" not in cons.columns:
            return float("nan")
        try:
            return float(cons.loc[(family, column, scope), "occupancy"])
        except (KeyError, TypeError):
            return float("nan")

    rows = []
    for r in calls.itertuples(index=False):
        rate_all = rate_for(r.family, r.column, SCOPE_ALL)
        rate_a = rate_for(r.family, r.column, r.group_a)
        rate_b = rate_for(r.family, r.column, r.group_b)
        occ = occupancy_for(r.family, r.column, SCOPE_ALL)
        rows.append(
            {
                "family": r.family,
                "partition": r.partition,
                "group_a": r.group_a,
                "group_b": r.group_b,
                "column": r.column,
                "type1": bool(r.type1),
                "type2": bool(r.type2),
                "rate_all": rate_all,
                "rate_a": rate_a,
                "rate_b": rate_b,
                "occupancy": occ,
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        out["class"] = pd.Series(dtype="string")
        out["status"] = pd.Series(dtype="string")
        return out

    missing = out[["rate_all", "rate_a", "rate_b"]].isna().any(axis=1)
    low_occupancy = out["occupancy"].notna() & (out["occupancy"] < t.min_occupancy)

    conserved_all = out["rate_all"] <= t.conserved_rate
    conserved_both = (out["rate_a"] <= t.conserved_rate) & (out["rate_b"] <= t.conserved_rate)
    variable_all = out["rate_all"] >= t.variable_rate
    divergent = out["type1"] | out["type2"]

    klass = np.select(
        [
            missing,
            conserved_all & divergent,
            out["type2"] & conserved_both,
            divergent & ~conserved_both,
            conserved_all & ~divergent,
            variable_all & ~divergent,
        ],
        [UNKNOWN, SUSPECT, DETERMINANT, RELAXED, CORE, VARIABLE],
        default=VARIABLE,
    )
    out["class"] = pd.Series(klass, index=out.index, dtype="string")

    status = np.select(
        [missing, low_occupancy],
        ["missing_conservation", "low_occupancy"],
        default="ok",
    )
    out["status"] = pd.Series(status, index=out.index, dtype="string")
    return out


def summary(classified: pd.DataFrame) -> pd.DataFrame:
    """Counts per class and status. Reported alongside results, never as a filter."""
    if classified.empty:
        return pd.DataFrame(columns=["family", "partition", "class", "status", "n"])
    return (
        classified.groupby(["family", "partition", "class", "status"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["family", "partition", "n"], ascending=[True, True, False])
        .reset_index(drop=True)
    )
