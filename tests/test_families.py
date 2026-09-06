"""Fusion and paralogue handling.

The fused case is modelled on a real one: a subtype whose adaptation module
declares two gene identities in a single slot, so one ORF carries both.
"""

import pandas as pd
import pytest

from contextshift.stages import families


def hits(rows):
    return pd.DataFrame(rows)


FUSED = hits([
    {"member_id": "orf1", "genome_id": "g1", "locus_id": "L1", "gene": "geneA", "start": 10},
    {"member_id": "orf1", "genome_id": "g1", "locus_id": "L1", "gene": "geneB", "start": 10},
    {"member_id": "orf2", "genome_id": "g1", "locus_id": "L1", "gene": "geneC", "start": 50},
])

TWO_COPIES = hits([
    {"member_id": "a", "genome_id": "g1", "locus_id": "L1", "gene": "geneB", "start": 10},
    {"member_id": "b", "genome_id": "g1", "locus_id": "L1", "gene": "geneB", "start": 90},
    {"member_id": "c", "genome_id": "g2", "locus_id": "L2", "gene": "geneB", "start": 10},
])


def test_fusion_is_flagged_and_kept_by_default():
    m, r = families.build(FUSED)
    assert r.n_fusion_orfs == 1
    assert r.fusion_genes == {"geneA+geneB": 1}
    assert r.n_excluded == 0
    # the fused ORF appears once per gene identity
    assert sorted(m[m["member_id"] == "orf1"]["family"]) == ["geneA", "geneB"]


def test_fusion_can_be_excluded_but_is_counted():
    m, r = families.build(FUSED, fusion_policy=families.EXCLUDE)
    assert r.n_excluded == 2
    assert r.n_fusion_orfs == 1
    assert "orf1" not in set(m["member_id"])
    assert set(m["member_id"]) == {"orf2"}


def test_unknown_policy_rejected():
    with pytest.raises(ValueError, match="fusion_policy"):
        families.build(FUSED, fusion_policy="drop-silently")


def test_two_copies_in_one_locus_are_not_collapsed():
    m, r = families.build(TWO_COPIES)
    assert len(m) == 3
    l1 = m[m["member_id"].isin(["a", "b"])].sort_values("copy_index")
    assert list(l1["copy_index"]) == [1, 2]
    assert r.n_multicopy_loci == 1
    assert r.max_copies == 2


def test_copy_index_follows_position_not_input_order():
    shuffled = TWO_COPIES.iloc[[1, 0, 2]].reset_index(drop=True)
    m, _ = families.build(shuffled)
    first = m[m["member_id"] == "a"]["copy_index"].iloc[0]
    second = m[m["member_id"] == "b"]["copy_index"].iloc[0]
    assert (first, second) == (1, 2)


def test_copies_in_separate_loci_both_index_one():
    m, _ = families.build(TWO_COPIES)
    assert int(m[m["member_id"] == "c"]["copy_index"].iloc[0]) == 1


def test_paralogues_under_distinct_labels_stay_separate():
    df = hits([
        {"member_id": "x", "genome_id": "g", "locus_id": "L", "gene": "cog1468", "start": 1},
        {"member_id": "y", "genome_id": "g", "locus_id": "L", "gene": "cog4343", "start": 5},
    ])
    m, r = families.build(df)
    assert set(r.families) == {"cog1468", "cog4343"}


def test_family_key_distinguishes_second_copy():
    m, _ = families.build(TWO_COPIES)
    keys = sorted(families.family_key(r) for r in m.itertuples(index=False))
    assert keys == ["geneB", "geneB", "geneB-2"]


def test_missing_columns_rejected():
    with pytest.raises(ValueError, match="missing columns"):
        families.build(pd.DataFrame({"member_id": ["a"]}))


def test_report_text_names_the_fused_combination():
    _, r = families.build(FUSED)
    assert "geneA+geneB" in r.as_text()
    assert "fusion ORFs     1" in r.as_text()


def test_split_keeps_unlabelled_members_visible():
    m, _ = families.build(TWO_COPIES)
    groups = families.split_by_partition(m, {"a": "one", "b": "two"})
    assert set(groups) == {"one", "two", "unlabelled"}
    assert list(groups["unlabelled"]["member_id"]) == ["c"]
