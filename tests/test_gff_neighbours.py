"""Gene neighbourhoods read from a GFF."""

import pytest

from contextshift.stages import neighbours

GFF = """##gff-version 3
c1\t.\tregion\t1\t9000\t.\t+\t.\tID=r1
c1\t.\tCDS\t100\t400\t.\t+\t0\tID=cds-p1;product=hypothetical protein
c1\t.\tCDS\t500\t800\t.\t+\t0\tID=cds-p2;product=widget synthase
c1\t.\tCDS\t900\t1200\t.\t-\t0\tID=cds-p3;product=target family protein
c1\t.\tCDS\t1300\t1600\t.\t+\t0\tID=cds-p4;product=marker protein Xyz1
c1\t.\tCDS\t1700\t2000\t.\t+\t0\tID=cds-p5;product=hypothetical protein
c2\t.\tCDS\t100\t400\t.\t+\t0\tID=cds-q1;product=target family protein
c2\t.\tCDS\t500\t800\t.\t+\t0\tID=cds-q2;product=hypothetical protein
"""


@pytest.fixture
def cds(tmp_path):
    p = tmp_path / "genomic.gff"
    p.write_text(GFF)
    return neighbours.read_gff_cds(p)


def test_reads_only_cds_ordered_along_each_contig(cds):
    assert list(cds["member_id"]) == ["p1", "p2", "p3", "p4", "p5", "q1", "q2"]
    assert list(cds[cds["contig"] == "c1"]["index_on_contig"]) == [0, 1, 2, 3, 4]
    assert list(cds[cds["contig"] == "c2"]["index_on_contig"]) == [0, 1]


def test_products_are_captured(cds):
    assert cds.set_index("member_id").loc["p4", "product"] == "marker protein Xyz1"


def test_neighbourhood_is_centred_and_signed(cds):
    near, truncated = neighbours.neighbourhood(cds, "p3", window=2)
    assert list(near["offset"]) == [-2, -1, 0, 1, 2]
    assert not truncated


def test_window_running_off_the_end_is_flagged(cds):
    _, truncated = neighbours.neighbourhood(cds, "p1", window=2)
    assert truncated
    _, truncated = neighbours.neighbourhood(cds, "q1", window=2)
    assert truncated


def test_neighbourhood_does_not_cross_contigs(cds):
    near, _ = neighbours.neighbourhood(cds, "q1", window=5)
    assert set(near["contig"]) == {"c2"}


def test_classification_finds_a_matching_neighbour(cds):
    out = neighbours.classify_by_neighbours(cds, ["p3"], pattern=r"marker", window=2)
    row = out.iloc[0]
    assert row["label"] == "associated"
    assert row["n_matching"] == 1
    assert int(row["nearest_offset"]) == 1


def test_a_member_with_no_match_is_solo(cds):
    out = neighbours.classify_by_neighbours(cds, ["p1"], pattern=r"marker", window=1)
    assert out.iloc[0]["label"] == "solo"


def test_edge_members_are_marked_so_a_negative_is_not_trusted(cds):
    out = neighbours.classify_by_neighbours(cds, ["q2"], pattern=r"marker", window=5)
    row = out.iloc[0]
    assert row["label"] == "solo"
    assert bool(row["edge"]), "a solo call at a contig edge must be flagged"


def test_the_member_itself_never_counts_as_its_own_neighbour(cds):
    out = neighbours.classify_by_neighbours(cds, ["p4"], pattern=r"marker", window=0)
    assert out.iloc[0]["label"] == "solo"


def test_a_member_absent_from_the_gff_is_reported_not_dropped(cds):
    out = neighbours.classify_by_neighbours(cds, ["nope"], pattern=r"marker")
    assert len(out) == 1
    assert out.iloc[0]["label"] == "not_in_gff"
