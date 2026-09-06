"""Per-site conservation, computed in-library or read from Rate4Site.

The default is the Jensen-Shannon divergence of each column against a
background amino-acid distribution, weighted by sequence redundancy
(Capra & Singh 2007, Bioinformatics 23:1875-1882), which that paper found to
outperform entropy and substitution-matrix scores at identifying functional
sites. It needs no external tool, so it is testable here.

Rate4Site (ConSurf's engine; Pupko et al. 2002) remains available where a
model-based evolutionary rate is wanted. Both are reported on the same scale
convention: LOWER means more constrained.

Run once per scope: the whole family, then each group alone.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

import pandas as pd

from ..schema import CONSERVATION
from ._external import RATE4SITE, run

SCOPE_ALL = "all"
_ROW = re.compile(r"^\s*(\d+)\s+(\S)\s+(-?[\d.eE+]+)")

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
GAPS = "-."

# BLOSUM62 background amino-acid frequencies, in AMINO_ACIDS order
BACKGROUND = (
    0.078, 0.019, 0.054, 0.063, 0.039, 0.074, 0.023, 0.052, 0.059, 0.096,
    0.024, 0.042, 0.052, 0.043, 0.051, 0.071, 0.059, 0.066, 0.014, 0.032,
)


def _sequence_weights(columns: list[str]) -> list[float]:
    """Henikoff position-based weights, normalised to mean 1."""
    n = len(columns[0]) if columns else 0
    weights = [0.0] * n
    for col in columns:
        seen = {c for c in col if c in AMINO_ACIDS}
        if not seen:
            continue
        counts = {a: col.count(a) for a in seen}
        for i, c in enumerate(col):
            if c in counts:
                weights[i] += 1.0 / (len(seen) * counts[c])
    total = sum(weights)
    return [w * n / total for w in weights] if total else [1.0] * n


def _jensen_shannon(freqs: list[float], background=BACKGROUND) -> float:
    import math

    m = [(p + q) / 2 for p, q in zip(freqs, background, strict=True)]
    def kl(a, b):
        return sum(x * math.log2(x / y) for x, y in zip(a, b, strict=True) if x > 0 and y > 0)
    return 0.5 * kl(freqs, m) + 0.5 * kl(background, m)


def jensen_shannon(
    alignment: Path,
    family: str,
    scope: str,
    pseudocount: float = 1e-6,
    gap_cutoff: float = 0.3,
) -> pd.DataFrame:
    """Per-column conservation as JS divergence from a background distribution.

    Returned as `rate` on the same convention as Rate4Site: lower is more
    constrained. Columns with more gaps than `gap_cutoff` get a NaN rate and an
    occupancy value, rather than a misleading score.
    """
    from Bio import AlignIO

    aln = AlignIO.read(str(alignment), "fasta")
    columns = [str(rec.seq).upper() for rec in aln]
    length = aln.get_alignment_length()
    weights = _sequence_weights([
        "".join(seq[i] for seq in columns) for i in range(length)
    ])

    rows = []
    for i in range(length):
        col = [seq[i] for seq in columns]
        occupancy = sum(1 for c in col if c not in GAPS) / len(col)
        if occupancy < (1 - gap_cutoff):
            rows.append({"column": i, "rate": float("nan"), "occupancy": occupancy})
            continue
        counts = [pseudocount] * len(AMINO_ACIDS)
        for c, w in zip(col, weights, strict=True):
            j = AMINO_ACIDS.find(c)
            if j >= 0:
                counts[j] += w
        total = sum(counts)
        freqs = [c / total for c in counts]
        divergence = _jensen_shannon(freqs)
        rows.append({"column": i, "rate": 1.0 - divergence, "occupancy": occupancy})

    df = pd.DataFrame(rows)
    df["family"] = family
    df["scope"] = scope
    return CONSERVATION.validate(df[["family", "column", "scope", "rate", "occupancy"]])


def rate4site(alignment: Path, family: str, scope: str, tree: Path | None = None) -> pd.DataFrame:
    RATE4SITE.require()
    tmp = Path(tempfile.mkdtemp(prefix="contextshift-r4s-"))
    out = tmp / "r4s.res"
    cmd = ["rate4site", "-s", str(alignment), "-o", str(out)]
    if tree is not None:
        cmd += ["-t", str(tree)]
    run(cmd, cwd=tmp)
    return parse_rate4site(out, family=family, scope=scope)


def parse_rate4site(path: Path, family: str, scope: str) -> pd.DataFrame:
    # layout unverified against real output; see VERIFICATION.md
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        m = _ROW.match(line)
        if m:
            rows.append({"column": int(m.group(1)), "rate": float(m.group(3))})
    df = pd.DataFrame(rows)
    df["family"] = family
    df["scope"] = scope
    return CONSERVATION.validate(df[["family", "column", "scope", "rate"]])


def with_occupancy(conservation: pd.DataFrame, occupancy: list[float]) -> pd.DataFrame:
    df = conservation.copy()
    lookup = {i + 1: v for i, v in enumerate(occupancy)}
    df["occupancy"] = df["column"].map(lookup).astype("float64")
    return CONSERVATION.validate(df)
