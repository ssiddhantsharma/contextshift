"""Column to residue mapping, checked against 4IC1 chain A coordinates."""

from pathlib import Path

import pytest

from contextshift.stages import mapping

CIF = Path(__file__).parent / "data" / "structures" / "4IC1_A.cif"
CHAIN = "A"

# SG atoms within 2.4 A of an Fe of the SF4 cluster, in all 10 chains of 4IC1
FE_S_CYSTEINES = [32, 188, 191, 197]
# heavy atoms within 2.3 A of the catalytic Mn, in all 10 chains
METAL_SITE = [62, 99, 113, 114]


@pytest.fixture(scope="module")
def residues():
    return mapping.read_chain(CIF, CHAIN)


def ungapped_alignment(residues):
    return residues.sequence


def test_author_numbering_may_be_negative(residues):
    # 4IC1 chain A is numbered from -3: an expression tag precedes residue 1.
    # Anything that assumes author numbering starts at 1, or equals a sequence
    # offset, is wrong here.
    assert residues.chain == CHAIN
    assert len(residues.resnums) == len(residues.resnames)
    assert residues.resnums == sorted(residues.resnums)
    assert min(residues.resnums) < 1
    assert 0 in residues.resnums


def test_round_trip_holds_on_a_gapless_alignment(residues):
    seq = ungapped_alignment(residues)
    m, mismatches = mapping.map_columns(seq, residues, "Cas4", "4IC1_A")
    assert not mismatches
    assert mapping.round_trip(m, seq, residues)


def test_mapping_recovers_the_iron_sulfur_cysteines(residues):
    seq = ungapped_alignment(residues)
    m, _ = mapping.map_columns(seq, residues, "Cas4", "4IC1_A")
    got = m[m["resnum"].isin(FE_S_CYSTEINES)]
    assert len(got) == len(FE_S_CYSTEINES)
    assert set(got["resname"]) == {"CYS"}


def test_mapping_recovers_the_metal_site(residues):
    seq = ungapped_alignment(residues)
    m, _ = mapping.map_columns(seq, residues, "Cas4", "4IC1_A")
    got = m[m["resnum"].isin(METAL_SITE)].sort_values("resnum")
    assert list(got["resname"]) == ["HIS", "ASP", "GLU", "ILE"]


def test_gaps_shift_columns_without_breaking_the_round_trip(residues):
    seq = ungapped_alignment(residues)
    gapped = "---" + seq[:50] + "-----" + seq[50:]
    sub = mapping.ChainResidues(CHAIN, residues.resnums, residues.resnames)
    m, mismatches = mapping.map_columns(gapped, sub, "Cas4", "4IC1_A")
    assert not mismatches
    assert mapping.round_trip(m, gapped, sub)
    # the first residue is no longer column 0
    assert int(m["column"].min()) == 3


def test_off_by_one_is_caught_not_tolerated(residues):
    # drop the first residue from the alignment: everything shifts by one
    seq = ungapped_alignment(residues)[1:]
    with pytest.raises(mapping.MappingError, match="ungapped residues"):
        mapping.map_columns(seq, residues, "Cas4", "4IC1_A")


def test_wrong_sequence_is_reported_not_silently_mapped(residues):
    seq = "A" * len(residues.resnums)
    with pytest.raises(mapping.MappingError, match="mismatches"):
        mapping.map_columns(seq, residues, "Cas4", "4IC1_A")
    _, mismatches = mapping.map_columns(seq, residues, "Cas4", "4IC1_A", strict=False)
    assert len(mismatches) > 100


def test_iron_sulfur_cysteines_are_the_closest_residues_to_the_cluster(residues):
    d = mapping.distance_to(CIF, CHAIN, residues, target_res_names=("SF4",))
    closest = set(d.nsmallest(4, "min_distance")["resnum"].astype(int))
    assert closest == set(FE_S_CYSTEINES)
    assert d[d["resnum"].isin(FE_S_CYSTEINES)]["min_distance"].max() < 2.5


def test_metal_site_residues_are_closest_to_the_manganese(residues):
    d = mapping.distance_to(CIF, CHAIN, residues, target_res_names=("MN",))
    closest = set(d.nsmallest(4, "min_distance")["resnum"].astype(int))
    assert closest == set(METAL_SITE)


def test_missing_target_raises_rather_than_returning_zeros(residues):
    with pytest.raises(mapping.MappingError, match="no target atoms"):
        mapping.distance_to(CIF, CHAIN, residues, target_res_names=("XYZ",))
