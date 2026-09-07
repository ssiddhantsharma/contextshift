"""Residue-to-ligand distances against published values.

The `.dist_to_lig` file is from the data supporting Capra & Singh 2007: a
header naming the PDB entry, chain and ligand, then one line per residue giving
the distance in angstroms to the nearest relevant ligand atom.

1A71 has no deposited hydrogens, so heavy-atom distances reproduce those values
exactly. Where hydrogens are deposited they do not: on 1ADB (1,410 H atoms) the
published distances run about one bond length shorter, and one residue is
published at 1.59 A, which is an H-to-H contact. This stage measures heavy-atom
distances deliberately.

Two things the comparison established about how to call `distance_to`. The
nearest ligand atom for a residue in one chain is often a copy of the ligand
bound to another, so restricting the target to a single chain gives the wrong
answer -- on 1ADB it moved the maximum deviation from 1.7 A to 33 A. The
fixture therefore keeps every copy of the ligand alongside chain A.
"""

import re
from pathlib import Path

import pytest

from contextshift.stages import mapping

DATA = Path(__file__).parent / "data" / "ligand"
CIF = DATA / "1A71_A.cif"
DIST = DATA / "1A71_A.dist_to_lig"
CHAIN = "A"


def published():
    lines = DIST.read_text().splitlines()
    header = re.match(r"#(\w{4})\s+(\S)\s+(\S+)", lines[0])
    assert header, lines[0]
    values = {}
    for line in lines[1:]:
        p = line.split()
        if len(p) == 3:
            values[int(p[0])] = (p[1], float(p[2]))
    return header.group(1).upper(), header.group(2), header.group(3), values


def test_header_names_the_entry_chain_and_ligand():
    pdb, chain, ligand, values = published()
    assert (pdb, chain, ligand) == ("1A71", "A", "NAD")
    assert len(values) > 300


def test_distances_reproduce_the_published_values_exactly():
    _, _, ligand, values = published()
    residues = mapping.read_chain(CIF, CHAIN)
    mine = mapping.distance_to(CIF, CHAIN, residues, target_res_names=(ligand,))

    compared = 0
    for row in mine.itertuples(index=False):
        num = int(row.resnum)
        if num not in values or row.min_distance != row.min_distance:
            continue
        compared += 1
        assert abs(values[num][1] - row.min_distance) < 0.01, (
            f"residue {num}: published {values[num][1]:.3f}, computed {row.min_distance:.3f}"
        )
    assert compared > 300


def test_residue_identity_agrees_with_the_published_file():
    """A distance is only meaningful if both sides mean the same residue."""
    _, _, _, values = published()
    residues = mapping.read_chain(CIF, CHAIN)
    letters = dict(zip(residues.resnums, residues.sequence, strict=True))
    checked = 0
    for num, (letter, _) in values.items():
        if num in letters:
            assert letters[num] == letter, f"residue {num}: {letters[num]} vs {letter}"
            checked += 1
    assert checked > 300


def test_a_missing_ligand_raises_rather_than_returning_zeros():
    residues = mapping.read_chain(CIF, CHAIN)
    with pytest.raises(mapping.MappingError, match="no target atoms"):
        mapping.distance_to(CIF, CHAIN, residues, target_res_names=("NOTALIGAND",))


def test_the_ligand_is_present_in_more_than_one_chain():
    """A residue's nearest ligand atom is often a copy bound to another chain."""
    from biotite.structure.io.pdbx import CIFFile, get_structure

    s = get_structure(CIFFile.read(str(CIF)), model=1)
    chains = sorted(set(s.chain_id[s.res_name == "NAD"]))
    assert len(chains) > 1, chains


def test_hydrogens_are_excluded():
    import biotite.structure as struc
    from biotite.structure.io.pdbx import CIFFile, get_structure

    s = get_structure(CIFFile.read(str(CIF)), model=1)
    assert int((s.element == "H").sum()) == 0, "fixture chosen for having no deposited H"
    assert struc.filter_amino_acids(s).sum() > 0
