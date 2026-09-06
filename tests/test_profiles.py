"""Profile-based family assignment.

Fixtures are the real Pfam Cas4 model (PF01930, 'Cas_Cas4' / DUF83) and twelve
RefSeq proteins whose titles say Cas4.
"""

from pathlib import Path

import pytest

from contextshift.stages import profiles

DATA = Path(__file__).parent / "data" / "profiles"
HMM = DATA / "PF01930.hmm"
FASTA = DATA / "cas4_refseq_sample.faa"


def test_fixture_is_the_expected_pfam_model():
    head = HMM.read_text().split("\n")[:3]
    assert head[0].startswith("HMMER3")
    assert "Cas_Cas4" in head[1]
    assert "PF01930" in head[2]


def test_scan_finds_cas4_in_cas4_sequences():
    hits = profiles.scan(FASTA, HMM)
    assert len(hits) > 0
    assert set(hits["profile"]) == {"Cas_Cas4"}
    assert (hits["evalue"] <= 1e-5).all()


def test_assignment_reports_margin_and_runner_up():
    assigned = profiles.assign(profiles.scan(FASTA, HMM), {"Cas_Cas4": "cas4"})
    assert set(assigned["family"]) == {"cas4"}
    # a single profile leaves no runner-up, so the margin is infinite
    assert all(m == float("inf") for m in assigned["margin"])
    assert set(assigned["runner_up"]) == {""}


def test_empty_hits_give_an_empty_assignment():
    import pandas as pd

    out = profiles.assign(pd.DataFrame(columns=["member_id", "profile", "score", "evalue"]))
    assert out.empty
    assert "margin" in out.columns


def test_ambiguous_flags_close_calls():
    import pandas as pd

    hits = pd.DataFrame([
        {"member_id": "a", "profile": "X", "profile_accession": "", "score": 100.0, "evalue": 1e-30},
        {"member_id": "a", "profile": "Y", "profile_accession": "", "score": 97.0, "evalue": 1e-29},
        {"member_id": "b", "profile": "X", "profile_accession": "", "score": 100.0, "evalue": 1e-30},
        {"member_id": "b", "profile": "Y", "profile_accession": "", "score": 20.0, "evalue": 1e-3},
    ])
    assigned = profiles.assign(hits)
    amb = profiles.ambiguous(assigned, min_margin=10.0)
    assert list(amb["member_id"]) == ["a"]


def test_a_title_is_not_evidence_of_family_membership(tmp_path):
    """RefSeq titles saying 'Csa1' are not Cas4-family Csa1.

    Of 150 bacterial RefSeq proteins titled 'Csa1 family protein', none hit
    PF01930 or PF06023 at E<=1e-3, and RefSeq holds zero archaeal ones -- for a
    family that is subtype I-A specific and largely archaeal. A profile hit is
    evidence; a title is not.
    """
    unrelated = tmp_path / "u.faa"
    unrelated.write_text(">plantlike\n" + "MEEALKRIAELESQVSSLTSQ" * 8 + "\n")
    hits = profiles.scan(unrelated, HMM)
    assert hits.empty


@pytest.mark.parametrize("evalue", [1e-10, 1e-5])
def test_threshold_is_applied(evalue):
    hits = profiles.scan(FASTA, HMM, evalue=evalue)
    assert (hits["evalue"] <= evalue).all()
