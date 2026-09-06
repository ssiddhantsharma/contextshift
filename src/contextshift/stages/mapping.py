"""Alignment column to structure residue, and distance to named features."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..schema import MAPPING

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}
GAPS = "-."


class MappingError(ValueError):
    pass


@dataclass(frozen=True)
class ChainResidues:
    chain: str
    resnums: list[int]
    resnames: list[str]

    @property
    def sequence(self) -> str:
        return "".join(THREE_TO_ONE.get(r, "X") for r in self.resnames)


def read_chain(path: Path, chain: str, model: int = 1) -> ChainResidues:
    """Observed protein residues of one chain, in author numbering."""
    import biotite.structure as struc
    from biotite.structure.io.pdbx import CIFFile, get_structure

    s = get_structure(CIFFile.read(str(path)), model=model)
    sel = s[(s.chain_id == chain) & struc.filter_amino_acids(s)]
    if len(sel) == 0:
        raise MappingError(f"chain {chain!r} has no amino acids in {path}")

    resnums, resnames = [], []
    for num, name in zip(sel.res_id, sel.res_name, strict=True):
        num = int(num)
        if not resnums or resnums[-1] != num:
            resnums.append(num)
            resnames.append(str(name))
    return ChainResidues(chain, resnums, resnames)


def read_alignment(path: Path) -> dict[str, str]:
    from Bio import AlignIO

    return {rec.id: str(rec.seq) for rec in AlignIO.read(str(path), "fasta")}


def map_columns(
    aligned_sequence: str,
    residues: ChainResidues,
    family: str,
    reference_id: str,
    strict: bool = True,
) -> tuple[pd.DataFrame, list[str]]:
    """Map each non-gap alignment column to an author residue number.

    Columns are 0-based, matching DIVERGE. Returns the mapping and any
    residue-letter mismatches; strict=True raises on mismatch.
    """
    ungapped = [(i, c) for i, c in enumerate(aligned_sequence) if c not in GAPS]
    if len(ungapped) != len(residues.resnums):
        raise MappingError(
            f"{reference_id}: {len(ungapped)} ungapped residues in the alignment but "
            f"{len(residues.resnums)} observed in chain {residues.chain}. "
            "Unobserved residues must be reconciled before mapping, not guessed."
        )

    rows, mismatches = [], []
    for (column, letter), num, name in zip(
        ungapped, residues.resnums, residues.resnames, strict=True
    ):
        expected = THREE_TO_ONE.get(name, "X")
        if letter.upper() != expected and expected != "X":
            mismatches.append(
                f"column {column}: alignment has {letter!r}, {residues.chain}{num} is {name}"
            )
        rows.append(
            {
                "family": family,
                "column": column,
                "reference_id": reference_id,
                "resnum": num,
                "resname": name,
            }
        )

    if mismatches and strict:
        raise MappingError(
            f"{len(mismatches)} alignment/structure mismatches; first three:\n  "
            + "\n  ".join(mismatches[:3])
        )
    return MAPPING.validate(pd.DataFrame(rows)), mismatches


def round_trip(mapping: pd.DataFrame, aligned_sequence: str, residues: ChainResidues) -> bool:
    """Column to resnum to residue letter must return the aligned letter."""
    by_num = dict(zip(residues.resnums, residues.resnames, strict=True))
    for r in mapping.itertuples(index=False):
        letter = aligned_sequence[int(r.column)]
        name = by_num[int(r.resnum)]
        if THREE_TO_ONE.get(name, "X") not in (letter.upper(), "X"):
            return False
    return True


def distance_to(
    path: Path,
    chain: str,
    residues: ChainResidues,
    target_chains: tuple[str, ...] = (),
    target_res_names: tuple[str, ...] = (),
    model: int = 1,
) -> pd.DataFrame:
    """Minimum heavy-atom distance from each residue to a target set.

    Targets are named by chain or by residue name. NaN where absent, never 0.
    """
    from biotite.structure.io.pdbx import CIFFile, get_structure

    s = get_structure(CIFFile.read(str(path)), model=model)
    s = s[s.element != "H"]

    mask = np.zeros(len(s), dtype=bool)
    if target_chains:
        mask |= np.isin(s.chain_id, list(target_chains))
    if target_res_names:
        mask |= np.isin(s.res_name, list(target_res_names))
    target = s[mask]
    if len(target) == 0:
        raise MappingError(
            f"no target atoms for chains={target_chains} res_names={target_res_names}"
        )

    source = s[s.chain_id == chain]
    rows = []
    for num in residues.resnums:
        atoms = source[source.res_id == num]
        if len(atoms) == 0:
            rows.append({"resnum": num, "min_distance": np.nan})
            continue
        d = np.linalg.norm(
            atoms.coord[:, None, :] - target.coord[None, :, :], axis=-1
        )
        rows.append({"resnum": num, "min_distance": float(d.min())})
    return pd.DataFrame(rows)


def annotate(sites: pd.DataFrame, mapping: pd.DataFrame, distances: pd.DataFrame) -> pd.DataFrame:
    """Attach resnum and distance to sites, keeping unmapped columns."""
    out = sites.merge(mapping[["family", "column", "resnum", "resname"]],
                      on=["family", "column"], how="left")
    return out.merge(distances, on="resnum", how="left")
