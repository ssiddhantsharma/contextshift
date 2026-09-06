"""Motif discovery, run per group rather than over the pooled family.

Run on a pooled alignable set, MEME mostly restates what the alignment already
shows. Run separately per group and compared, it can find a motif private to
one group; that comparison is what this stage is for. Width is not fixed,
because real motif widths are not.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

import pandas as pd

from ..schema import MOTIFS
from ._external import MEME, run

_MOTIF = re.compile(
    r"^MOTIF\s+(\S+)\s+MEME-(\d+)\s+width\s*=\s*(\d+)\s+sites\s*=\s*(\d+)\s+llr\s*=\s*\S+\s+E-value\s*=\s*(\S+)",
    re.MULTILINE,
)


def meme(
    fasta: Path,
    family: str,
    partition: str,
    group: str,
    nmotifs: int = 10,
    minw: int = 6,
    maxw: int = 30,
    mod: str = "zoops",
) -> pd.DataFrame:
    MEME.require()
    tmp = Path(tempfile.mkdtemp(prefix="contextshift-meme-"))
    run(
        [
            "meme", str(fasta), "-protein", "-oc", str(tmp),
            "-nmotifs", str(nmotifs), "-minw", str(minw), "-maxw", str(maxw), f"-{mod}",
        ]
    )
    return parse_meme(tmp / "meme.txt", family, partition, group)


def parse_meme(path: Path, family: str, partition: str, group: str) -> pd.DataFrame:
    text = Path(path).read_text()
    rows = [
        {
            "family": family,
            "partition": partition,
            "group": group,
            "motif_id": m.group(2),
            "consensus": m.group(1),
            "width": int(m.group(3)),
            "n_sites": int(m.group(4)),
            "evalue": float(m.group(5)),
        }
        for m in _MOTIF.finditer(text)
    ]
    df = pd.DataFrame(rows, columns=[c.name for c in MOTIFS.columns])
    return MOTIFS.validate(df)


def private_motifs(motifs: pd.DataFrame) -> pd.DataFrame:
    """Motifs whose consensus occurs in exactly one group of a partition."""
    df = MOTIFS.validate(motifs.copy())
    spread = df.groupby(["family", "partition", "consensus"])["group"].nunique()
    spread.name = "n_groups"
    merged = df.merge(spread.reset_index(), on=["family", "partition", "consensus"])
    return merged[merged["n_groups"] == 1].reset_index(drop=True)
