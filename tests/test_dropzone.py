import pandas as pd

from contextshift.reference import Claim, DropZone, Slot, evaluate
from contextshift.reference.dropzone import Comparison
from contextshift.schema import MEMBERS


def compare_counts(computed, reference) -> Comparison:
    overlap = len(set(computed) & set(reference["member_id"]))
    return Comparison(
        slot="members",
        metric="overlap",
        value=overlap / max(len(reference), 1),
        expected=1.0,
        passed=overlap > 0,
    )


def make_zone(tmp_path):
    zone = DropZone(root=tmp_path)
    zone.register(
        Slot(
            name="members",
            filename="members.parquet",
            description="member table from the prior study",
            schema=MEMBERS,
            comparator=compare_counts,
        )
    )
    zone.register(Slot(name="motifs", filename="motifs.txt", description="motif definitions"))
    return zone


def test_absent_slots_are_reported_not_skipped(tmp_path):
    zone = make_zone(tmp_path)
    assert len(zone.absent()) == 2
    assert "[ ] members" in zone.report()
    assert "motif definitions" in zone.report()


def test_duplicate_slot_rejected(tmp_path):
    zone = make_zone(tmp_path)
    try:
        zone.register(Slot(name="members", filename="x", description="dup"))
    except ValueError as exc:
        assert "already registered" in str(exc)
    else:
        raise AssertionError("duplicate slot was accepted")


def test_present_slot_compares_without_replacing(tmp_path):
    zone = make_zone(tmp_path)
    MEMBERS.write(
        pd.DataFrame({"member_id": ["a", "b"], "genome_id": ["g", "g"], "family": ["F", "F"]}),
        tmp_path / "members.parquet",
    )
    computed = ["a", "b", "c"]
    comparisons = zone.compare({"members": computed})
    assert len(comparisons) == 1
    assert comparisons[0].value == 1.0
    assert computed == ["a", "b", "c"]


def test_unmapped_files_are_surfaced(tmp_path):
    zone = make_zone(tmp_path)
    (tmp_path / "surprise.csv").write_text("x\n")
    assert [p.name for p in zone.unmapped()] == ["surprise.csv"]
    assert "surprise.csv" in zone.report()


def test_refuted_claim_records_observation_rather_than_raising():
    claim = Claim(
        name="motif is group specific",
        statement="motif occurs in one group only",
        source="prior study",
        check=lambda ctx: (ctx["n_groups"] == 1, ctx["n_groups"]),
    )
    [result] = evaluate([claim], {"n_groups": 4})
    assert result.verdict == "REFUTED"
    assert result.observed == 4


def test_untestable_claim_is_an_error_not_a_pass():
    claim = Claim(
        name="broken",
        statement="",
        source="",
        check=lambda ctx: (_ for _ in ()).throw(KeyError("missing")),
    )
    [result] = evaluate([claim], {})
    assert result.verdict == "ERROR"
    assert not result.upheld
