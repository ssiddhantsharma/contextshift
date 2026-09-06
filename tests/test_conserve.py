"""In-library conservation (Capra & Singh 2007 Jensen-Shannon)."""

import math

import pytest

from contextshift.stages import conserve


def write_alignment(tmp_path, seqs):
    p = tmp_path / "a.fa"
    p.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)))
    return p


def test_invariant_column_is_more_constrained_than_a_variable_one(tmp_path):
    p = write_alignment(tmp_path, ["AK", "AR", "AD", "AE", "AW", "AY"])
    d = conserve.jensen_shannon(p, "F", "all").set_index("column")
    assert d.loc[0, "rate"] < d.loc[1, "rate"]


def test_rate_convention_is_lower_means_conserved(tmp_path):
    p = write_alignment(tmp_path, ["CCCC", "CCCA", "CCCD", "CCCE"])
    d = conserve.jensen_shannon(p, "F", "all").set_index("column")
    assert d.loc[0, "rate"] == pytest.approx(d.loc[1, "rate"])
    assert d.loc[0, "rate"] < d.loc[3, "rate"]


def test_gappy_columns_get_nan_not_a_misleading_score(tmp_path):
    p = write_alignment(tmp_path, ["A-", "A-", "A-", "AK"])
    d = conserve.jensen_shannon(p, "F", "all").set_index("column")
    assert math.isnan(d.loc[1, "rate"])
    assert d.loc[1, "occupancy"] == pytest.approx(0.25)
    assert not math.isnan(d.loc[0, "rate"])


def test_occupancy_is_reported_for_every_column(tmp_path):
    p = write_alignment(tmp_path, ["AK", "A-", "AK", "AK"])
    d = conserve.jensen_shannon(p, "F", "all")
    assert d["occupancy"].notna().all()


def test_redundant_sequences_do_not_dominate(tmp_path):
    # nine identical sequences plus one different: weighting should stop the
    # duplicated residue from being scored as near-invariant
    d1 = tmp_path / "one"; d1.mkdir()
    d2 = tmp_path / "two"; d2.mkdir()
    balanced = conserve.jensen_shannon(write_alignment(d1, ["AK", "AD"]), "F", "all")
    skewed = conserve.jensen_shannon(write_alignment(d2, ["AK"] * 9 + ["AD"]), "F", "all")
    b = balanced.set_index("column").loc[1, "rate"]
    s = skewed.set_index("column").loc[1, "rate"]
    assert abs(b - s) < 0.25, (b, s)


def test_scope_and_family_are_carried(tmp_path):
    p = write_alignment(tmp_path, ["AK", "AR"])
    d = conserve.jensen_shannon(p, "Fam", "group_x")
    assert set(d["family"]) == {"Fam"}
    assert set(d["scope"]) == {"group_x"}


def test_rate4site_parser_still_available():
    assert callable(conserve.parse_rate4site)


def test_columns_are_zero_based(tmp_path):
    p = write_alignment(tmp_path, ["AK", "AR"])
    d = conserve.jensen_shannon(p, "F", "all")
    assert int(d["column"].min()) == 0
    assert int(d["column"].max()) == 1
