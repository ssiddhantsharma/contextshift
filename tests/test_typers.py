"""Parsing three typers, and comparing their partitions."""

import pytest

from contextshift.partition import adjusted_rand_index
from contextshift.stages import loci

DF_HEADER = "sys_id\ttype\tsubtype\tsys_beg\tsys_end\tprotein_in_syst\tgenes_count\tname_of_profiles_in_sys"


def df_file(tmp_path, rows):
    p = tmp_path / "defense_finder_systems.tsv"
    p.write_text(DF_HEADER + "\n" + "\n".join(rows) + "\n")
    return p


def test_defensefinder_subtype_string_is_parsed(tmp_path):
    p = df_file(tmp_path, [
        "c1_CAS_1\tCAS\tCAS_Class1-Subtype-I-E\tp1\tp3\tp1,p2,p3\t3\tCas1,Cas2,Cas3",
        "c2_CAS_1\tCAS\tCAS_Class1-Subtype-I-C\tp4\tp5\tp4,p5\t2\tCas1,Cas4",
    ])
    out = loci.parse_defensefinder(p, system="CAS")
    assert set(out["subtype"]) == {"I-E", "I-C"}
    assert set(out["member_id"]) == {"p1", "p2", "p3", "p4", "p5"}


def test_profiles_are_never_zipped_to_proteins(tmp_path):
    """Both lists are alphabetised independently; pairing them mislabels genes."""
    p = df_file(tmp_path, [
        "c1_CAS_1\tCAS\tCAS_Class1-Subtype-I-A\tp1\tp3\tzz,aa,mm\t3\tCas4,Cas1,Cas2",
    ])
    out = loci.parse_defensefinder(p)
    assert "gene" not in out.columns
    assert set(out["member_id"]) == {"zz", "aa", "mm"}


def test_a_system_pattern_selects_one_family(tmp_path):
    p = df_file(tmp_path, [
        "c1_RM_1\tRM\tRM_type_I\tp1\tp2\tp1,p2\t2\tRM_1,RM_2",
        "c1_CAS_1\tCAS\tCAS_Class1-Subtype-I-B\tp3\tp4\tp3,p4\t2\tCas1,Cas4",
    ])
    assert set(loci.parse_defensefinder(p, system="CAS")["subtype"]) == {"I-B"}
    assert len(loci.parse_defensefinder(p)) == 4


def test_defensefinder_missing_columns_rejected(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("sys_id\ttype\n a\tb\n")
    with pytest.raises(loci.LocusError, match="missing columns"):
        loci.parse_defensefinder(p)


def padloc_file(tmp_path, rows):
    p = tmp_path / "x_padloc.csv"
    header = "system.number,seqid,system,target.name,hmm.accession,hmm.name,protein.name,start,end,strand"
    p.write_text(header + "\n" + "\n".join(rows) + "\n")
    return p


def test_padloc_dotted_columns_and_subtype(tmp_path):
    p = padloc_file(tmp_path, [
        "1,c1,cas_type_I-E,p1,PF001,h1,Cas1,1,900,+",
        "1,c1,cas_type_I-E,p2,PF002,h2,Cas2,901,1200,+",
    ])
    out = loci.parse_padloc(p)
    assert list(out["member_id"]) == ["p1", "p2"]
    assert set(out["subtype"]) == {"I-E"}
    assert set(out["locus_id"]) == {"c1@1"}


def test_padloc_system_pattern(tmp_path):
    p = padloc_file(tmp_path, [
        "1,c1,RM_type_I,p9,PF9,h9,RM,1,50,+",
        "2,c1,cas_type_I-C,p1,PF1,h1,Cas4,1,900,+",
    ])
    assert set(loci.parse_padloc(p, system="cas")["subtype"]) == {"I-C"}


def test_partitions_from_two_typers_are_comparable(tmp_path):
    dfp = df_file(tmp_path, [
        "c1_CAS_1\tCAS\tCAS_Class1-Subtype-I-E\tp1\tp2\tp1,p2\t2\tCas1,Cas2",
        "c2_CAS_1\tCAS\tCAS_Class1-Subtype-I-C\tp3\tp4\tp3,p4\t2\tCas1,Cas4",
    ])
    pad = padloc_file(tmp_path, [
        "1,c1,cas_type_I-E,p1,PF1,h1,Cas1,1,9,+",
        "1,c1,cas_type_I-E,p2,PF2,h2,Cas2,10,20,+",
        "2,c2,cas_type_I-C,p3,PF3,h3,Cas1,1,9,+",
        "2,c2,cas_type_I-C,p4,PF4,h4,Cas4,10,20,+",
    ])
    a = loci.partition_from(loci.parse_defensefinder(dfp), "subtype", "defensefinder-3", "df")
    b = loci.partition_from(loci.parse_padloc(pad), "subtype", "padloc-2", "padloc")
    assert adjusted_rand_index(a, b) == pytest.approx(1.0)


def test_disagreement_shows_up_as_a_lower_index(tmp_path):
    dfp = df_file(tmp_path, [
        "c1_CAS_1\tCAS\tCAS_Class1-Subtype-I-E\tp1\tp4\tp1,p2,p3,p4\t4\tCas1,Cas2,Cas3,Cas5",
    ])
    pad = padloc_file(tmp_path, [
        "1,c1,cas_type_I-E,p1,PF1,h1,Cas1,1,9,+",
        "1,c1,cas_type_I-E,p2,PF2,h2,Cas2,10,20,+",
        "2,c1,cas_type_I-C,p3,PF3,h3,Cas1,21,30,+",
        "2,c1,cas_type_I-C,p4,PF4,h4,Cas4,31,40,+",
    ])
    a = loci.partition_from(loci.parse_defensefinder(dfp), "subtype", "df", "s")
    b = loci.partition_from(loci.parse_padloc(pad), "subtype", "padloc", "s")
    assert adjusted_rand_index(a, b) < 1.0


def test_padloc_subdivided_subtypes_are_preserved(tmp_path):
    """PADLOC splits I-B into I-B1/I-B2; that is real, not a parse error."""
    p = padloc_file(tmp_path, [
        "1,c1,cas_type_I-B1,p1,PF1,h1,Cas4,1,9,+",
        "2,c2,cas_type_I-B2,p2,PF2,h2,Cas4,1,9,+",
    ])
    out = loci.parse_padloc(p)
    assert set(out["subtype"]) == {"I-B1", "I-B2"}


def test_coarse_subtype_only_for_cross_scheme_comparison():
    assert loci.coarse_subtype("I-B1") == "I-B"
    assert loci.coarse_subtype("I-F3") == "I-F"
    assert loci.coarse_subtype("I-E") == "I-E"
    assert loci.coarse_subtype("CAS_Class1-Subtype-I-E") == "CAS_Class1-Subtype-I-E"


def test_padloc_adaptation_system_has_no_subtype(tmp_path):
    """cas_adaptation is a real PADLOC system with no subtype; keep it visible."""
    p = padloc_file(tmp_path, [
        "1,c1,cas_adaptation,p1,PF1,h1,Cas4,1,9,+",
    ])
    out = loci.parse_padloc(p)
    assert out["subtype"].iloc[0] == "cas_adaptation"


def test_granularity_difference_is_reconcilable(tmp_path):
    """Coarsening should raise agreement when the only difference is depth."""
    dfp = df_file(tmp_path, [
        "c1_CAS_1\tCAS\tCAS_Class1-Subtype-I-B\tp1\tp2\tp1,p2\t2\tCas1,Cas4",
    ])
    pad = padloc_file(tmp_path, [
        "1,c1,cas_type_I-B1,p1,PF1,h1,Cas1,1,9,+",
        "1,c1,cas_type_I-B1,p2,PF2,h2,Cas4,10,20,+",
    ])
    a = loci.partition_from(loci.parse_defensefinder(dfp), "subtype", "df", "s")
    raw = loci.parse_padloc(pad)
    coarse = raw.assign(subtype=raw["subtype"].map(loci.coarse_subtype))
    b = loci.partition_from(coarse, "subtype", "padloc", "s")
    assert a.labels == b.labels
