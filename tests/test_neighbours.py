"""Neighbourhood parsing and group comparison."""

import pytest

from contextshift.partition import Partition, Provenance
from contextshift.stages import neighbours

PROV = Provenance(source="test")


def operon_row(query, species, cluster, neighbour, start):
    return "\t".join([
        f"{query}|{species}", "300", "+", "+", str(cluster),
        "0", "100", str(start), str(start + 100), f"{neighbour}#1",
        "GCF_1", "assembly", "info",
    ])


@pytest.fixture
def operon_file(tmp_path):
    rows = []
    # two members that share clusters 5 and 6 around themselves
    for q in ("Q1", "Q2"):
        rows.append(operon_row(q, "Sp", 5, "N5", 100))
        rows.append(operon_row(q, "Sp", 0, q, 200))
        rows.append(operon_row(q, "Sp", 6, "N6", 300))
    # two members whose neighbours are all singletons
    for q, n in (("Q3", "X3"), ("Q4", "X4")):
        rows.append(operon_row(q, "Sp", 0, n, 100))
        rows.append(operon_row(q, "Sp", 0, q, 200))
    p = tmp_path / "x_operon.tsv"
    p.write_text("\n".join(rows) + "\n")
    return p


@pytest.fixture
def partition():
    return Partition(
        name="context",
        labels={"Q1": "in_system", "Q2": "in_system", "Q3": "orphan", "Q4": "orphan"},
        provenance=PROV,
    )


def test_parses_positional_columns(operon_file):
    n = neighbours.parse_operon_tsv(operon_file)
    assert set(n["member_id"]) == {"Q1", "Q2", "Q3", "Q4"}
    assert set(n["neighbour_id"]) >= {"N5", "N6", "X3", "X4"}


def test_query_sits_at_offset_zero(operon_file):
    n = neighbours.parse_operon_tsv(operon_file)
    q1 = n[n["member_id"] == "Q1"].sort_values("offset")
    self_row = q1[q1["neighbour_id"] == "Q1"]
    assert int(self_row["offset"].iloc[0]) == 0
    assert list(q1["offset"]) == [-1, 0, 1]


def test_conserved_clusters_found_for_the_shared_group(operon_file, partition):
    n = neighbours.parse_operon_tsv(operon_file)
    c = neighbours.conservation(n, partition)
    kept = c[(c["group"] == "in_system") & c["conserved"]]
    assert set(kept["cluster_id"]) == {"5", "6"}


def test_group_with_no_shared_context_reports_none(operon_file, partition):
    n = neighbours.parse_operon_tsv(operon_file)
    summary = neighbours.shared_and_private(neighbours.conservation(n, partition))
    assert summary["n_conserved"]["orphan"] == 0
    assert "orphan" in summary["groups_with_none"]


def test_unclustered_neighbours_are_not_counted_as_conserved(operon_file, partition):
    n = neighbours.parse_operon_tsv(operon_file)
    c = neighbours.conservation(n, partition)
    assert not c[(c["cluster_id"] == "0") & c["conserved"]].shape[0]


def test_private_clusters_are_attributed_to_one_group(operon_file, partition):
    n = neighbours.parse_operon_tsv(operon_file)
    summary = neighbours.shared_and_private(neighbours.conservation(n, partition))
    assert set(summary["private"]["in_system"]) == {"5", "6"}
    assert summary["shared"] == []


def test_members_outside_the_partition_are_excluded_not_mislabelled(operon_file):
    n = neighbours.parse_operon_tsv(operon_file)
    small = Partition(name="context", labels={"Q1": "a"}, provenance=PROV)
    c = neighbours.conservation(n, small)
    assert set(c["group"]) == {"a"}


def test_annotation_join(operon_file, tmp_path):
    desc = tmp_path / "d.txt"
    desc.write_text("1(2)\tN5\tsome conserved protein\n")
    n = neighbours.parse_operon_tsv(operon_file)
    out = neighbours.annotate(n, neighbours.parse_outdesc(desc))
    assert out[out["neighbour_id"] == "N5"]["annotation"].iloc[0] == "some conserved protein"
    assert out[out["neighbour_id"] == "N6"]["annotation"].isna().all()


def test_outdesc_parser_reads_the_shipped_layout(tmp_path):
    p = tmp_path / "o.txt"
    p.write_text("1(2)\tWP_049155312.1\tMULTISPECIES: glycogen debranching protein GlgX \n")
    assert neighbours.parse_outdesc(p) == {
        "WP_049155312.1": "MULTISPECIES: glycogen debranching protein GlgX"
    }
