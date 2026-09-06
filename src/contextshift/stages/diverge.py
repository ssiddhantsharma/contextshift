"""Adapter for DIVERGE v4 (Gu et al., MBE 2025; zjupgx/diverge4).

DIVERGE ships a Python API and a Streamlit UI but no CLI, and it is strict
about tree shape: Newick, no internal node names, depth >= 3, one tree per
cluster. Meeting those constraints reproducibly is most of the work, so tree
conformance is done here and tested independently of whether diverge4 is
installed.

The result-parsing layer has NOT yet been validated against an installed
diverge4; `normalise` is written against the documented column meanings and
must be checked against real output before any result is trusted.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..partition import Partition
from ..schema import SITES

TYPE1 = "type1"
TYPE2 = "type2"

_INTERNAL_LABEL = re.compile(r"\)[^(),:;]+(?=[:,);])")


class DivergeUnavailable(RuntimeError):
    pass


@dataclass(frozen=True)
class TreeCheck:
    ok: bool
    depth: int
    n_leaves: int
    problems: tuple[str, ...] = ()


def strip_internal_labels(newick: str) -> str:
    """Remove internal node names. DIVERGE rejects trees that carry them."""
    return _INTERNAL_LABEL.sub(")", newick)


def _depth(newick: str) -> int:
    depth = best = 0
    for ch in newick:
        if ch == "(":
            depth += 1
            best = max(best, depth)
        elif ch == ")":
            depth -= 1
    return best


def check_tree(newick: str, min_depth: int = 3, min_leaves: int = 4) -> TreeCheck:
    problems: list[str] = []
    depth = _depth(newick)
    n_leaves = newick.count(",") + 1 if newick.strip() else 0
    if depth < min_depth:
        problems.append(f"depth {depth} < {min_depth}")
    if n_leaves < min_leaves:
        problems.append(f"{n_leaves} leaves < {min_leaves}")
    if _INTERNAL_LABEL.search(newick):
        problems.append("internal node labels present")
    return TreeCheck(not problems, depth, n_leaves, tuple(problems))


def conform(newick: str, min_depth: int = 3, min_leaves: int = 4) -> tuple[str, TreeCheck]:
    cleaned = strip_internal_labels(newick.strip())
    if not cleaned.endswith(";"):
        cleaned += ";"
    return cleaned, check_tree(cleaned, min_depth, min_leaves)


def write_cluster_trees(
    trees: dict[str, str], outdir: Path, min_depth: int = 3, min_leaves: int = 4
) -> tuple[dict[str, Path], dict[str, TreeCheck]]:
    """Write one conformant Newick per group. Rejects are returned, not dropped."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    checks: dict[str, TreeCheck] = {}
    for group, newick in trees.items():
        cleaned, check = conform(newick, min_depth, min_leaves)
        checks[group] = check
        if check.ok:
            path = outdir / f"{group.replace('/', '_')}.nwk"
            path.write_text(cleaned + "\n")
            paths[group] = path
    return paths, checks


def _import_diverge():
    try:
        import diverge
    except ImportError as exc:
        raise DivergeUnavailable(
            "diverge4 is not installed; `uv pip install diverge4` "
            "(see github.com/zjupgx/diverge4)"
        ) from exc
    return diverge


def available() -> bool:
    try:
        _import_diverge()
    except DivergeUnavailable:
        return False
    return True


def normalise(
    raw: pd.DataFrame,
    family: str,
    partition: str,
    group_a: str,
    group_b: str,
    test: str,
) -> pd.DataFrame:
    """Map a DIVERGE result frame onto the SITES schema.

    UNVERIFIED against installed diverge4 output. Column names below are the
    documented ones; confirm before relying on any number this produces.
    """
    df = raw.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]

    def pick(*names: str) -> pd.Series:
        for n in names:
            if n in df.columns:
                return df[n]
        return pd.Series([pd.NA] * len(df), index=df.index, dtype="float64")

    out = pd.DataFrame(
        {
            "family": family,
            "partition": partition,
            "group_a": group_a,
            "group_b": group_b,
            "column": pick("position", "site", "column").astype("Int64"),
            "test": test,
            "theta": pd.to_numeric(pick("theta", "theta_ml"), errors="coerce"),
            "statistic": pd.to_numeric(pick("lrt", "statistic", "z"), errors="coerce"),
            "posterior": pd.to_numeric(pick("posterior", "qk"), errors="coerce"),
            "pvalue": pd.to_numeric(pick("pvalue", "p_value", "p"), errors="coerce"),
        }
    )
    return SITES.validate(out)


def run_pair(
    alignment: Path,
    tree_a: Path,
    tree_b: Path,
    family: str,
    partition: str,
    group_a: str,
    group_b: str,
) -> pd.DataFrame:
    """Run Gu99 (Type-I) and the Type-II test for one pair of groups."""
    diverge = _import_diverge()
    frames = []

    gu = diverge.Gu99(str(alignment), trees=[str(tree_a), str(tree_b)])
    frames.append(normalise(gu.results(), family, partition, group_a, group_b, TYPE1))

    type2 = getattr(diverge, "Type2", None)
    if type2 is not None:
        t2 = type2(str(alignment), trees=[str(tree_a), str(tree_b)])
        frames.append(normalise(t2.results(), family, partition, group_a, group_b, TYPE2))

    return pd.concat(frames, ignore_index=True)


def run_partition(
    alignment: Path,
    trees: dict[str, Path],
    partition: Partition,
    family: str,
    skip_underpowered: bool = False,
) -> tuple[pd.DataFrame, list[str]]:
    """Run every group pair. Skipped pairs are returned, never silently dropped."""
    skipped: list[str] = []
    weak = set(partition.underpowered_groups())
    frames = []
    for a, b in partition.pairs():
        if a not in trees or b not in trees:
            skipped.append(f"{a}|{b}: missing conformant tree")
            continue
        if skip_underpowered and (a in weak or b in weak):
            skipped.append(f"{a}|{b}: underpowered")
            continue
        frames.append(run_pair(alignment, trees[a], trees[b], family, partition.name, a, b))

    if not frames:
        return SITES.validate(pd.DataFrame(columns=[c.name for c in SITES.columns])), skipped
    return pd.concat(frames, ignore_index=True), skipped
