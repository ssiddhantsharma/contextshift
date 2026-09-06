"""Typed-locus parsing. List columns arrive as stringified Python lists."""

import pandas as pd
import pytest

from contextshift.stages import loci

HEADER = "Contig\tOperon\tStart\tEnd\tPrediction\tGenes\tPositions\tE-values"


def row(contig, operon, pred, genes, positions, start=1):
    return "\t".join([
        contig, operon, str(start), str(start + 900), pred,
        str(genes), str(positions), str(["1e-20"] * len(genes)),
    ])


@pytest.fixture
def operons_file(tmp_path):
    p = tmp_path / "cas_operons.tab"
    p.write_text("\n".join([
        HEADER,
        row("c1", "c1@1", "I-A", ["Cas1_", "Cas2", "Cas4"], ["p1", "p2", "p3"]),
        row("c2", "c2@1", "I-E", ["Cas1_", "Cas2"], ["p4", "p5"]),
        row("c3", "c3@1", "Ambiguous", ["Cas1_"], ["p6"]),
    ]) + "\n")
    return p


def test_list_columns_are_parsed_not_left_as_strings(operons_file):
    df = loci.read_operons(operons_file)
    assert isinstance(df["Genes"].iloc[0], list)
    assert df["Genes"].iloc[0] == ["Cas1_", "Cas2", "Cas4"]
    assert df["Positions"].iloc[0] == ["p1", "p2", "p3"]


def test_missing_columns_rejected(tmp_path):
    p = tmp_path / "bad.tab"
    p.write_text("Contig\tOperon\n c\to\n")
    with pytest.raises(loci.LocusError, match="missing columns"):
        loci.read_operons(p)


def test_hits_are_one_row_per_gene(operons_file):
    hits = loci.to_hits(loci.read_operons(operons_file))
    assert len(hits) == 6
    assert list(hits[hits["locus_id"] == "c1@1"]["gene"]) == ["Cas1_", "Cas2", "Cas4"]
    assert list(hits[hits["locus_id"] == "c1@1"]["member_id"]) == ["p1", "p2", "p3"]


def test_parallel_lists_of_unequal_length_are_rejected():
    df = pd.DataFrame([{
        "Contig": "c", "Operon": "o", "Start": 1, "End": 2, "Prediction": "I-A",
        "Genes": ["a", "b"], "Positions": ["p1"],
    }])
    with pytest.raises(loci.LocusError, match="2 genes but 1 positions"):
        loci.to_hits(df)


def test_subtype_partition_records_the_scheme(operons_file):
    hits = loci.to_hits(loci.read_operons(operons_file))
    p, dropped = loci.subtype_partition(hits, "typer-1.8.0", "some-scheme-2020")
    assert p.name == "subtype"
    assert p.provenance.scheme == "some-scheme-2020"
    assert set(p.group_names) == {"I-A", "I-E"}


def test_ambiguous_members_are_returned_not_silently_dropped(operons_file):
    hits = loci.to_hits(loci.read_operons(operons_file))
    _, dropped = loci.subtype_partition(hits, "t", "s")
    assert dropped == ["p6"]


def test_ambiguous_can_be_kept(operons_file):
    hits = loci.to_hits(loci.read_operons(operons_file))
    p, dropped = loci.subtype_partition(hits, "t", "s", drop_ambiguous=False)
    assert "Ambiguous" in p.group_names
    assert dropped == []


def test_all_ambiguous_raises_rather_than_returning_empty(tmp_path):
    p = tmp_path / "a.tab"
    p.write_text("\n".join([HEADER, row("c", "c@1", "False", ["Cas1_"], ["p1"])]) + "\n")
    hits = loci.to_hits(loci.read_operons(p))
    with pytest.raises(loci.LocusError, match="no members left"):
        loci.subtype_partition(hits, "t", "s")


def test_gene_presence_is_measured_per_subtype(operons_file):
    df = loci.read_operons(operons_file)
    pres = loci.gene_presence(df, "Cas4")
    got = dict(zip(pres["subtype"], pres["fraction"], strict=True))
    assert got["I-A"] == 1.0
    assert got["I-E"] == 0.0


def test_hits_feed_the_families_stage(operons_file):
    from contextshift.stages import families

    hits = loci.to_hits(loci.read_operons(operons_file))
    members, report = families.build(hits)
    assert report.n_members == 6
    assert set(report.families) == {"Cas1_", "Cas2", "Cas4"}
