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


@pytest.mark.parametrize(
    "seqs", [["WWWW", "WWWW"], ["ACDE", "FGHI"], ["AAAA", "YYYY"]]
)
def test_jensen_shannon_is_bounded(tmp_path, seqs):
    """With log base 2 the divergence lies in [0, 1] (Capra & Singh 2007)."""
    d = conserve.jensen_shannon(write_alignment(tmp_path, seqs), "F", "all")
    assert d["rate"].min() >= 0.0
    assert d["rate"].max() <= 1.0


def test_identical_distributions_have_zero_divergence():
    background = list(conserve.BACKGROUND)
    assert conserve._jensen_shannon(background, background) == pytest.approx(0.0, abs=1e-12)


def test_divergence_is_symmetric():
    skewed = [0.5, 0.5] + [0.0] * 18
    background = list(conserve.BACKGROUND)
    assert conserve._jensen_shannon(skewed, background) == pytest.approx(
        conserve._jensen_shannon(background, skewed)
    )


def test_a_background_column_is_less_constrained_than_an_invariant_one(tmp_path):
    counts = [max(1, round(f * 200)) for f in conserve.BACKGROUND]
    column = "".join(a * n for a, n in zip(conserve.AMINO_ACIDS, counts, strict=True))

    mixed = tmp_path / "mixed"; mixed.mkdir()
    single = tmp_path / "single"; single.mkdir()
    background = conserve.jensen_shannon(write_alignment(mixed, list(column)), "F", "all")
    invariant = conserve.jensen_shannon(
        write_alignment(single, ["W"] * len(column)), "F", "all"
    )
    assert invariant.loc[0, "rate"] < background.loc[0, "rate"]
