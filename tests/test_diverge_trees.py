"""DIVERGE rejects trees on shape, so conformance is checked without it installed."""

from contextshift.stages.diverge import (
    check_tree,
    conform,
    strip_internal_labels,
    write_cluster_trees,
)

DEEP = "(((a:0.1,b:0.1):0.2,(c:0.1,d:0.1):0.2):0.3,e:0.4);"
LABELLED = "(((a:0.1,b:0.1)node1:0.2,(c:0.1,d:0.1)node2:0.2)root:0.3,e:0.4);"
FLAT = "(a:0.1,b:0.1,c:0.1);"


def test_strips_internal_labels():
    assert "node1" not in strip_internal_labels(LABELLED)
    assert "node2" not in strip_internal_labels(LABELLED)


def test_strip_keeps_leaf_names():
    out = strip_internal_labels(LABELLED)
    for leaf in "abcde":
        assert f"{leaf}:" in out


def test_labelled_tree_fails_check_before_conforming():
    assert not check_tree(LABELLED).ok


def test_conform_makes_it_pass():
    cleaned, check = conform(LABELLED)
    assert check.ok
    assert cleaned.endswith(";")


def test_shallow_tree_rejected():
    check = check_tree(FLAT)
    assert not check.ok
    assert any("depth" in p for p in check.problems)


def test_depth_measured():
    assert check_tree(DEEP).depth == 3


def test_write_returns_rejects_rather_than_dropping(tmp_path):
    paths, checks = write_cluster_trees({"good": DEEP, "bad": FLAT}, tmp_path)
    assert set(paths) == {"good"}
    assert set(checks) == {"good", "bad"}
    assert not checks["bad"].ok
