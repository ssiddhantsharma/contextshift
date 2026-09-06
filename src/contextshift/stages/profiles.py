"""Assign family membership by profile rather than by a typer's gene label.

A gene name from a typer, or a database record's title, is a convenience; a
profile hit is evidence. This matters most where one loose label spans proteins
that two profiles separate, and where an unrelated family happens to share the
label. Titles have been observed to do both.

HMMs come from InterPro's API, which serves them gzipped. pfam.xfam.org is
retired and returns nothing.
"""

from __future__ import annotations

import gzip
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

INTERPRO_HMM = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{acc}/?annotation=hmm"


class ProfileError(ValueError):
    pass


@dataclass(frozen=True)
class Profile:
    accession: str
    label: str


def fetch_pfam(accession: str, outdir: Path, timeout: int = 60) -> Path:
    """Download one Pfam HMM. Cached by accession."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{accession}.hmm"
    if path.exists() and path.stat().st_size:
        return path

    url = INTERPRO_HMM.format(acc=accession)
    with urllib.request.urlopen(url, timeout=timeout) as r:
        raw = r.read()
    text = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
    if not text.startswith(b"HMMER3"):
        raise ProfileError(f"{accession}: response is not an HMMER3 profile")
    path.write_bytes(text)
    return path


def combine(paths: list[Path], out: Path) -> Path:
    out = Path(out)
    out.write_bytes(b"".join(Path(p).read_bytes() for p in paths))
    return out


def scan(fasta: Path, hmm: Path, cpus: int = 0, evalue: float = 1e-5) -> pd.DataFrame:
    """hmmsearch every profile against every sequence."""
    import pyhmmer
    from pyhmmer import easel, plan7

    with easel.SequenceFile(str(fasta), digital=True) as sf:
        sequences = sf.read_block()

    def text(v) -> str:
        # pyhmmer returned bytes before 0.11 and str after
        return v.decode() if isinstance(v, bytes) else ("" if v is None else str(v))

    rows = []
    with plan7.HMMFile(str(hmm)) as hf:
        for hits in pyhmmer.hmmsearch(hf, sequences, cpus=cpus, E=evalue):
            profile = text(hits.query.name)
            accession = text(hits.query.accession)
            for hit in hits:
                if hit.evalue > evalue:
                    continue
                rows.append({
                    "member_id": text(hit.name),
                    "profile": profile,
                    "profile_accession": accession,
                    "score": float(hit.score),
                    "evalue": float(hit.evalue),
                })
    return pd.DataFrame(
        rows, columns=["member_id", "profile", "profile_accession", "score", "evalue"]
    )


def assign(hits: pd.DataFrame, profiles: dict[str, str] | None = None) -> pd.DataFrame:
    """Best-scoring profile per sequence, with the runner-up kept.

    The margin between best and second is reported so a near-tie between two
    profiles is visible rather than resolved silently.
    """
    if hits.empty:
        return pd.DataFrame(
            columns=["member_id", "family", "profile", "score", "evalue", "margin", "runner_up"]
        )
    ordered = hits.sort_values(["member_id", "score"], ascending=[True, False])
    rows = []
    for member, grp in ordered.groupby("member_id", sort=False):
        best = grp.iloc[0]
        second = grp.iloc[1] if len(grp) > 1 else None
        label = (profiles or {}).get(best["profile"], best["profile"])
        rows.append({
            "member_id": member,
            "family": label,
            "profile": best["profile"],
            "score": float(best["score"]),
            "evalue": float(best["evalue"]),
            "margin": float(best["score"] - second["score"]) if second is not None else float("inf"),
            "runner_up": second["profile"] if second is not None else "",
        })
    return pd.DataFrame(rows)


def ambiguous(assigned: pd.DataFrame, min_margin: float = 10.0) -> pd.DataFrame:
    """Sequences whose top two profiles are too close to separate."""
    return assigned[assigned["margin"] < min_margin].reset_index(drop=True)
