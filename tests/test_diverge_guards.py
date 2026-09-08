"""Guards against three measured DIVERGE 4.1.0 defects.

The boundaries here are measured, not guessed. Substituting one internal branch
length in the package's own cl1.tree: 0.0 runs, 1e-7 through 1e-5 abort the
interpreter with SIGABRT, 1e-4 and above run.
"""

import pandas as pd
import pytest

from contextshift.stages.diverge import (
    MIN_SAFE_BRANCH,
    TYPE1,
    TYPE2,
    Convergence,
    assess,
    conform,
    floor_branches,
    gap_free_columns,
    unsafe_branches,
)

SAFE = "(((a:0.1,b:0.1):0.1,c:0.1):0.1,d:0.1);"
TINY = "(((a:0.1,b:0.000002):0.1,c:0.1):0.1,d:0.1);"
ZERO = "(((a:0.1,b:0.0):0.1,c:0.1):0.1,d:0.1);"


def test_aborting_range_is_detected():
    assert unsafe_branches(TINY)[0] == 1
    assert unsafe_branches(TINY)[1] == pytest.approx(2e-6)


def test_exactly_zero_is_left_alone():
    # 0.0 is accepted by DIVERGE; only (0, 1e-4) aborts
    assert unsafe_branches(ZERO)[0] == 0
    assert floor_branches(ZERO)[1] == 0


def test_safe_lengths_untouched():
    assert unsafe_branches(SAFE)[0] == 0
    assert floor_branches(SAFE) == (SAFE, 0)


@pytest.mark.parametrize("length,unsafe", [
    (0.0, 0), (1e-7, 1), (1e-6, 1), (1e-5, 1), (1e-4, 0), (1e-3, 0),
])
def test_measured_boundary(length, unsafe):
    tree = f"(((a:0.1,b:{length:.10f}):0.1,c:0.1):0.1,d:0.1);"
    assert unsafe_branches(tree)[0] == unsafe


def test_flooring_is_idempotent():
    once, n = floor_branches(TINY)
    assert n == 1
    assert unsafe_branches(once)[0] == 0
    assert floor_branches(once)[1] == 0


def test_conform_floors_and_records():
    cleaned, check = conform(TINY)
    assert check.floored == 1
    assert unsafe_branches(cleaned)[0] == 0


def test_conform_can_leave_branches_alone():
    cleaned, check = conform(TINY, floor=None)
    assert check.floored == 0
    assert unsafe_branches(cleaned)[0] == 1


def test_gap_free_columns(tmp_path):
    fasta = tmp_path / "a.aln"
    fasta.write_text(">x\nAC-DE\n>y\nACG-E\n>z\nACGDE\n")
    # column 0,1,4 are gap-free; 2 and 3 each carry a gap
    assert gap_free_columns(fasta) == 3


def _frames(mfe, ml, se, lrt, posteriors):
    comparisons = pd.DataFrame(
        [{"test": TYPE1, "parameter": k, "value": v}
         for k, v in (("MFE Theta", mfe), ("ThetaML", ml),
                      ("SE Theta", se), ("LRT Theta", lrt))]
        + [{"test": TYPE2, "parameter": "Alpha ML", "value": 1.0}]
    )
    sites = pd.DataFrame({"test": [TYPE1] * len(posteriors), "posterior": posteriors})
    return comparisons, sites


def test_healthy_matches_the_reference_signature():
    # the CASP fixture: 0.156 / 0.124 with SE and LRT present and posteriors set
    verdict = assess(*_frames(0.156, 0.124, 0.030, 16.8, [0.1, 0.2, 0.8]))
    assert verdict.verdict == "healthy"
    assert verdict.ok


def test_collapsed_ml_estimate_fails():
    verdict = assess(*_frames(0.407, 0.0, None, None, [None, None]))
    assert verdict.verdict == "failed"
    assert not verdict.ok
    assert any("collapsed" in p for p in verdict.problems)


def test_all_nan_posteriors_fail():
    verdict = assess(*_frames(0.3, 0.3, 0.05, 10.0, [None, None, None]))
    assert verdict.verdict == "failed"
    assert any("NaN" in p for p in verdict.problems)


def test_missing_variance_step_is_degraded_not_failed():
    # Cas4's signature: estimators agree, posteriors present, SE and LRT absent
    verdict = assess(*_frames(0.352, 0.352, None, None, [0.1, 0.9]))
    assert verdict.verdict == "degraded"
    assert verdict.ok


def test_disagreeing_estimators_fail():
    verdict = assess(*_frames(0.60, 0.10, 0.05, 10.0, [0.1, 0.2]))
    assert verdict.verdict == "failed"
    assert any("disagree" in p for p in verdict.problems)


def test_min_safe_branch_is_the_measured_value():
    assert MIN_SAFE_BRANCH == 1e-4


def test_convergence_ok_property():
    assert Convergence("healthy").ok
    assert Convergence("degraded", ("x",)).ok
    assert not Convergence("failed", ("x",)).ok
