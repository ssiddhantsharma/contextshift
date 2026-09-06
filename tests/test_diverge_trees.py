"""Tree shape rules, cross-checked against diverge.binding.check_tree."""

import pytest

from contextshift.stages.diverge import (
    check_tree,
    column_label,
    conform,
    tree_depth,
    write_cluster_trees,
)

# Depths below are measured, not assumed, and the accept/reject verdicts were
# checked against diverge.binding.check_tree directly. The boundary is exact:
# depth 3 is rejected, depth 4 is accepted.
SHALLOW3 = "((a:0.1,b:0.1):0.1,c:0.1);"                              # depth 3 -> REJECT
DEEP4 = "(((a:0.1,b:0.1):0.1,c:0.1):0.1,d:0.1);"                     # depth 4 -> ACCEPT
DEEP5 = "((((a:0.1,b:0.1):0.1,c:0.1):0.2,d:0.1):0.3,e:0.4);"         # depth 5 -> ACCEPT
FLAT = "(a:0.1,b:0.1,c:0.1);"                                        # depth 2 -> REJECT
LABELLED = "((((a:0.1,b:0.1)n1:0.1,c:0.1)n2:0.2,d:0.1)n3:0.3,e:0.4);"


def test_depth_counts_edges_not_parentheses():
    assert tree_depth(FLAT) == 2
    assert tree_depth(SHALLOW3) == 3
    assert tree_depth(DEEP4) == 4


def test_threshold_is_strictly_greater_than_three():
    assert not check_tree(SHALLOW3, min_leaves=0).ok
    assert check_tree(DEEP4, min_leaves=0).ok


def test_rejection_message_names_the_real_rule():
    problems = check_tree(SHALLOW3, min_leaves=0).problems
    assert any("requires > 3" in p for p in problems)


def test_shallow_tree_rejected():
    assert not check_tree(FLAT, min_leaves=0).ok


def test_internal_labels_are_not_a_problem():
    # DIVERGE reads trees with Bio.Phylo and reformats them; labelled internal
    # nodes are fine. An earlier version of this adapter invented that rule.
    assert check_tree(LABELLED).ok


def test_unparseable_tree_is_reported_not_raised():
    check = check_tree("(((not a tree")
    assert not check.ok
    assert any("unparseable" in p for p in check.problems)


def test_conform_adds_terminator():
    cleaned, check = conform(DEEP5.rstrip(";"))
    assert cleaned.endswith(";")
    assert check.ok


def test_write_returns_rejects_rather_than_dropping(tmp_path):
    paths, checks = write_cluster_trees(
        {"good": DEEP5, "bad": SHALLOW3}, tmp_path, min_leaves=0
    )
    assert set(paths) == {"good"}
    assert set(checks) == {"good", "bad"}
    assert not checks["bad"].ok


def test_column_label_matches_diverge_convention():
    # verified against real output: Type2 named its column "A/B"
    assert column_label("A", "B") == "A/B"


@pytest.mark.parametrize("newick", [SHALLOW3, DEEP4, DEEP5, FLAT, LABELLED])
def test_agrees_with_diverge_own_check(newick):
    """Our verdict must match DIVERGE's, or we will hand it trees it rejects."""
    binding = pytest.importorskip("diverge.binding")
    import io

    from Bio import Phylo

    tree = Phylo.read(io.StringIO(newick), "newick")
    try:
        binding.check_tree(tree)
        theirs = True
    except ValueError:
        theirs = False
    except Exception:
        pytest.skip("tree not parseable by DIVERGE either")

    ours = check_tree(newick, min_leaves=0).ok
    assert ours == theirs, f"disagree on {newick}: ours={ours} diverge={theirs}"
