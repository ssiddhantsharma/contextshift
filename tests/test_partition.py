import math

import pytest

from contextshift.partition import Partition, Provenance, adjusted_rand_index

PROV = Provenance(source="test", scheme="none")


def make(labels, weights=None):
    return Partition(name="p", labels=labels, provenance=PROV, weights=weights or {})


def test_rejects_empty():
    with pytest.raises(ValueError):
        make({})


def test_rejects_weights_for_unknown_members():
    with pytest.raises(ValueError):
        make({"a": "x"}, weights={"b": 1.0})


def test_counts_and_groups():
    p = make({"a": "x", "b": "x", "c": "y"})
    assert p.counts() == {"x": 2, "y": 1}
    assert p.groups() == {"x": ["a", "b"], "y": ["c"]}
    assert p.pairs() == [("x", "y")]


def test_effective_n_equals_n_when_unweighted():
    p = make({f"m{i}": "x" for i in range(10)})
    assert p.effective_counts()["x"] == pytest.approx(10.0)


def test_effective_n_collapses_redundant_members():
    # one genuinely independent sequence plus nine near-identical copies of another
    weights = {"m0": 1.0} | {f"m{i}": 1 / 9 for i in range(1, 10)}
    p = make({f"m{i}": "x" for i in range(10)}, weights=weights)
    eff = p.effective_counts()["x"]
    assert eff < 10.0
    assert eff == pytest.approx(2.0, abs=0.5)


def test_underpowered_groups_flagged():
    p = make({f"m{i}": "small" for i in range(5)})
    assert p.underpowered_groups() == ["small"]


def test_cross_makes_cells():
    a = make({"m1": "solo", "m2": "assoc"})
    b = Partition(
        name="q", labels={"m1": "I-A", "m2": "I-C"}, provenance=PROV
    )
    c = a.cross(b)
    assert c.labels == {"m1": "solo|I-A", "m2": "assoc|I-C"}
    assert c.name == "p|q"


def test_cross_without_shared_members_raises():
    a = make({"m1": "x"})
    b = Partition(name="q", labels={"m9": "y"}, provenance=PROV)
    with pytest.raises(ValueError):
        a.cross(b)


def test_restrict_and_drop():
    p = make({"a": "x", "b": "y", "c": "y"})
    assert p.restrict(["a", "b"]).counts() == {"x": 1, "y": 1}
    assert p.drop_groups(["y"]).counts() == {"x": 1}


def test_ari_identical_is_one():
    a = make({"a": "1", "b": "1", "c": "2", "d": "2"})
    b = make({"a": "x", "b": "x", "c": "y", "d": "y"})
    assert adjusted_rand_index(a, b) == pytest.approx(1.0)


def test_ari_uncorrelated_is_near_zero():
    a = make({f"m{i}": str(i % 2) for i in range(40)})
    b = make({f"m{i}": str((i // 20) % 2) for i in range(40)})
    assert abs(adjusted_rand_index(a, b)) < 0.2


def test_ari_on_disjoint_members_is_nan():
    a = make({"a": "x"})
    b = make({"z": "y"})
    assert math.isnan(adjusted_rand_index(a, b))


def test_json_roundtrip(tmp_path):
    p = make({"a": "x", "b": "y"}, weights={"a": 0.5})
    path = tmp_path / "p.json"
    p.to_json(path)
    back = Partition.from_json(path)
    assert back.labels == p.labels
    assert back.weights == p.weights
    assert back.provenance == p.provenance


def test_tsv_with_weights(tmp_path):
    path = tmp_path / "p.tsv"
    path.write_text("# comment\na\tsolo\t0.25\nb\tassoc\n")
    p = Partition.from_tsv(path, name="context", provenance=PROV)
    assert p.labels == {"a": "solo", "b": "assoc"}
    assert p.weight("a") == 0.25
    assert p.weight("b") == 1.0


def test_caveats_travel_with_power():
    p = Partition(
        name="p",
        labels={f"m{i}": ("clean" if i < 40 else "messy") for i in range(80)},
        provenance=PROV,
        caveats={"messy": "labelling does not correspond to a clade"},
    )
    rows = {r.label: r for r in p.power()}
    assert rows["messy"].caveat
    assert rows["messy"].qualified
    assert not rows["clean"].qualified


def test_caveat_for_unknown_group_rejected():
    with pytest.raises(ValueError, match="unknown groups"):
        Partition(name="p", labels={"a": "x"}, provenance=PROV, caveats={"nope": "hi"})


def test_caveats_survive_restrict_and_roundtrip(tmp_path):
    p = Partition(
        name="p", labels={"a": "x", "b": "y"}, provenance=PROV,
        caveats={"x": "polyphyletic"},
    )
    assert p.restrict(["a"]).caveats == {"x": "polyphyletic"}
    path = tmp_path / "p.json"
    p.to_json(path)
    assert Partition.from_json(path).caveats == {"x": "polyphyletic"}
