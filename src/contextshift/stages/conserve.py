"""Per-site evolutionary rate via Rate4Site, ConSurf's engine.

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
