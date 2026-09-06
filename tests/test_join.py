"""The 2x2 is the result, so its corners are pinned here."""

import pandas as pd

from contextshift.join import CORE, DETERMINANT, RELAXED, SUSPECT, UNKNOWN, classify, summary


def site(column, test, posterior):
    """DIVERGE reports a posterior Qk in [0, 1]; higher means more divergent."""
    return {
        "family": "F",
        "partition": "context",
        "group_a": "A",
        "group_b": "B",
        "column": column,
        "test": test,
        "posterior": posterior,
    }


def cons(column, scope, rate, occupancy=1.0):
    return {"family": "F", "column": column, "scope": scope, "rate": rate, "occupancy": occupancy}


def build(sites, conservation):
    return classify(pd.DataFrame(sites), pd.DataFrame(conservation))


def test_catalytic_core_is_not_called_divergent():
    # conserved everywhere, no significant test: the negative control
    out = build(
        [site(10, "type1", 0.10), site(10, "type2", 0.05)],
        [cons(10, "all", 0.05), cons(10, "A", 0.05), cons(10, "B", 0.05)],
    )
    assert out.loc[0, "class"] == CORE


def test_specificity_determinant():
    # constrained inside each group, property shift between them
    out = build(
        [site(20, "type2", 0.97)],
        [cons(20, "all", 1.2), cons(20, "A", 0.1), cons(20, "B", 0.1)],
    )
    assert out.loc[0, "class"] == DETERMINANT


def test_relaxed_constraint_is_type1_not_determinant():
    out = build(
        [site(30, "type1", 0.97)],
        [cons(30, "all", 1.0), cons(30, "A", 0.1), cons(30, "B", 2.0)],
    )
    assert out.loc[0, "class"] == RELAXED


def test_conserved_everywhere_but_divergent_is_suspect():
    out = build(
        [site(40, "type2", 0.97)],
        [cons(40, "all", 0.05), cons(40, "A", 0.05), cons(40, "B", 0.05)],
    )
    assert out.loc[0, "class"] == SUSPECT


def test_missing_conservation_is_surfaced_not_dropped():
    out = build([site(50, "type2", 0.97)], [cons(50, "all", 0.2)])
    assert len(out) == 1
    assert out.loc[0, "class"] == UNKNOWN
    assert out.loc[0, "status"] == "missing_conservation"


def test_low_occupancy_is_flagged_not_filtered():
    out = build(
        [site(60, "type2", 0.97)],
        [cons(60, "all", 1.2, occupancy=0.1), cons(60, "A", 0.1), cons(60, "B", 0.1)],
    )
    assert len(out) == 1
    assert out.loc[0, "status"] == "low_occupancy"
    assert out.loc[0, "class"] == DETERMINANT


def test_posterior_below_threshold_is_not_a_determinant():
    sites = pd.DataFrame([site(70, "type2", 0.40)])
    out = classify(
        sites,
        pd.DataFrame([cons(70, "all", 1.2), cons(70, "A", 0.1), cons(70, "B", 0.1)]),
    )
    assert out.loc[0, "class"] != DETERMINANT


def test_every_tested_column_appears_once():
    sites = [site(c, "type2", 0.20) for c in range(1, 6)]
    conservation = [
        row for c in range(1, 6) for row in (cons(c, "all", 1.0), cons(c, "A", 1.0), cons(c, "B", 1.0))
    ]
    out = build(sites, conservation)
    assert len(out) == 5
    assert sorted(out["column"]) == [1, 2, 3, 4, 5]


def test_summary_counts_classes():
    out = build(
        [site(10, "type2", 0.05), site(20, "type2", 0.97)],
        [
            cons(10, "all", 0.05), cons(10, "A", 0.05), cons(10, "B", 0.05),
            cons(20, "all", 1.2), cons(20, "A", 0.1), cons(20, "B", 0.1),
        ],
    )
    s = summary(out)
    assert set(s["class"]) == {CORE, DETERMINANT}
    assert s["n"].sum() == 2


def test_empty_input_returns_empty_frame():
    out = classify(
        pd.DataFrame(columns=["family", "partition", "group_a", "group_b", "column", "test",
                              "posterior"]),
        pd.DataFrame(columns=["family", "column", "scope", "rate"]),
    )
    assert out.empty
