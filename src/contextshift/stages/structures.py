"""Predicted structures from the AlphaFold database.

Solved structures are scarce for many families, so mapping determinants across
groups needs predictions. Two things are checked before a model is used: the
sequence it was built from must match the sequence being analysed, and its
confidence must clear a floor. A model that fails either is reported, not
silently mapped.

The file URL carries a model version that changes (v3 and v4 now 404; v6 is
current), so it is read from the API rather than constructed.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

#: Mean pLDDT below which backbone geometry should not carry an argument.
MIN_PLDDT = 70.0


class StructureError(ValueError):
    pass


@dataclass(frozen=True)
class Model:
    accession: str
    version: int
    plddt: float
    sequence: str
    cif_url: str
    path: Path | None = None

    @property
    def confident(self) -> bool:
        return self.plddt >= MIN_PLDDT


def lookup(accession: str, timeout: int = 45) -> Model | None:
    """Metadata for one UniProt accession, or None if AFDB has no model."""
    try:
        with urllib.request.urlopen(API.format(acc=accession), timeout=timeout) as r:
            payload = json.load(r)
    except urllib.error.HTTPError as exc:
        # 404 is "no model"; 400 is a malformed accession. Neither is fatal to
        # a batch, so both return None and are counted by the caller.
        if exc.code in (400, 404):
            return None
        raise
    if not payload:
        return None
    d = payload[0]
    return Model(
        accession=d.get("uniprotAccession", accession),
        version=int(d.get("latestVersion", 0)),
        plddt=float(d.get("globalMetricValue", 0.0)),
        sequence=d.get("uniprotSequence", ""),
        cif_url=d.get("cifUrl", ""),
    )


def fetch(accession: str, outdir: Path, timeout: int = 120) -> Model | None:
    """Download the mmCIF for one accession. Cached by accession and version."""
    model = lookup(accession, timeout=timeout)
    if model is None:
        return None
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"AF-{model.accession}-v{model.version}.cif"
    if not (path.exists() and path.stat().st_size):
        with urllib.request.urlopen(model.cif_url, timeout=timeout) as r:
            path.write_bytes(r.read())
    return Model(**{**model.__dict__, "path": path})


def check(model: Model, expected_sequence: str, min_plddt: float = MIN_PLDDT) -> list[str]:
    """Reasons this model should not be used. Empty means it is usable."""
    problems = []
    if model.plddt < min_plddt:
        problems.append(f"mean pLDDT {model.plddt:.1f} < {min_plddt}")
    if expected_sequence and model.sequence != expected_sequence:
        if len(model.sequence) != len(expected_sequence):
            problems.append(
                f"sequence length {len(model.sequence)} != {len(expected_sequence)}"
            )
        else:
            n = sum(a != b for a, b in zip(model.sequence, expected_sequence, strict=True))
            problems.append(f"sequence differs at {n} positions")
    return problems


def fetch_many(
    accessions: list[str], outdir: Path, expected: dict[str, str] | None = None
) -> tuple[pd.DataFrame, list[str]]:
    """Fetch models for many accessions. Every miss and rejection is returned."""
    rows, notes = [], []
    for acc in accessions:
        try:
            model = fetch(acc, outdir)
        except Exception as exc:
            notes.append(f"{acc}: fetch failed ({exc})")
            continue
        if model is None:
            notes.append(f"{acc}: no AlphaFold model")
            continue
        problems = check(model, (expected or {}).get(acc, ""))
        for p in problems:
            notes.append(f"{acc}: {p}")
        rows.append({
            "accession": model.accession,
            "version": model.version,
            "plddt": model.plddt,
            "path": str(model.path),
            "usable": not problems,
        })
    return pd.DataFrame(
        rows, columns=["accession", "version", "plddt", "path", "usable"]
    ), notes
