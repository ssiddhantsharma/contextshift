"""Adapter for DIVERGE v4 (Gu et al., MBE 2025; zjupgx/diverge4).

Verified against diverge 4.1.0 built from source and run on the package's own
CASP test data, not against its documentation. Four things the docs get wrong
or leave implicit, each of which silently corrupts results:

1. Tree files are POSITIONAL arguments. The `trees=` keyword takes Bio.Phylo
   Tree objects, not paths, so passing paths there fails or misbehaves.
2. `.summary` and `.results` are PROPERTIES, not methods.
3. `.results` carries the alignment position in its INDEX (named "Position"),
   not in a column, and the index is sparse: only positions DIVERGE kept
   survive. On the CASP data 781 of 2088 columns came back, indexed 132..1460.
4. The value is a POSTERIOR PROBABILITY Qk in [0, 1], not a p-value. Theta and
   alpha are per-comparison and live in `.summary`, not per site.

Upstream defect in 4.1.0: `Gu99`, `Rvs` and `TypeOneAnalysis` call an
undefined `get_colnames` and raise NameError on construction, which makes every
Type-I entry point unusable as shipped. The calculators already return
display-ready labels, so the missing function is a pass-through; `apply_shim`
installs it and Gu99 then runs, returning 781 kept positions on CASP with the
same index as Type2. The shim is opt-out and always recorded, never silent.
`Gu2001` failed independently on this build and is not used.
"""

from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..partition import Partition
from ..schema import COMPARISONS, SITES

TYPE1 = "type1"
TYPE2 = "type2"

POSITION_BASE = 0

DEFAULT_QK = 0.9

BROKEN_UPSTREAM = {
    "Gu99": "calls undefined get_colnames() in .summary (diverge 4.1.0)",
    "Rvs": "calls undefined get_colnames() in .summary (diverge 4.1.0)",
    "TypeOneAnalysis": "calls undefined get_colnames() in .summary (diverge 4.1.0)",
}

SHIM_NOTE = (
    "Type-I ran through a shim for diverge 4.1.0's undefined get_colnames "
    "(a pass-through over calculator._r_names())"
)


class DivergeUnavailable(RuntimeError):
    pass


class DivergeUpstreamBug(RuntimeError):
    pass


@dataclass(frozen=True)
class TreeCheck:
    ok: bool
    depth: int
    n_leaves: int
    problems: tuple[str, ...] = ()


def _read_tree(newick: str):
    from Bio import Phylo

    return Phylo.read(io.StringIO(newick), "newick")


def tree_depth(newick: str) -> int:
    """Edges from root to deepest leaf, as DIVERGE measures it."""
    tree = _read_tree(newick)
    return max(len(tree.trace(tree.root, clade)) for clade in tree.get_terminals())


def check_tree(newick: str, min_depth: int = 3, min_leaves: int = 4) -> TreeCheck:
    """DIVERGE requires depth strictly greater than min_depth."""
    problems: list[str] = []
    try:
        tree = _read_tree(newick)
        depth = max(len(tree.trace(tree.root, c)) for c in tree.get_terminals())
        n_leaves = len(tree.get_terminals())
    except Exception as exc:
        return TreeCheck(False, 0, 0, (f"unparseable: {exc}",))

    if depth <= min_depth:
        problems.append(f"depth {depth} <= {min_depth} (DIVERGE requires > {min_depth})")
    if n_leaves < min_leaves:
        problems.append(f"{n_leaves} leaves < {min_leaves}")
    return TreeCheck(not problems, depth, n_leaves, tuple(problems))


def conform(newick: str, min_depth: int = 3, min_leaves: int = 4) -> tuple[str, TreeCheck]:
    cleaned = newick.strip()
    if not cleaned.endswith(";"):
        cleaned += ";"
    return cleaned, check_tree(cleaned, min_depth, min_leaves)


def write_cluster_trees(
    trees: dict[str, str], outdir: Path, min_depth: int = 3, min_leaves: int = 4
) -> tuple[dict[str, Path], dict[str, TreeCheck]]:
    """Write one Newick per group; rejects are returned, not dropped."""
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
            "diverge is not installed. The PyPI name is `diverge`, not `diverge4`, "
            "its sdist is broken (setup.py reads a requirements.txt it does not ship), "
            "and the default src/ tree is Windows-targeted. Build from a clone with "
            "src_linux/ instead. See github.com/zjupgx/diverge4"
        ) from exc
    return diverge


def needs_shim() -> bool:
    """True if this build lacks get_colnames (diverge 4.1.0 defect)."""
    try:
        import diverge.binding as binding
    except ImportError:
        return False
    return not hasattr(binding, "get_colnames")


def apply_shim() -> bool:
    """Install the missing get_colnames. True if it was needed."""
    if not needs_shim():
        return False
    import diverge.binding as binding

    binding.get_colnames = lambda names: list(names)
    return True


def available() -> bool:
    try:
        _import_diverge()
    except DivergeUnavailable:
        return False
    return True


def column_label(group_a: str, group_b: str) -> str:
    """DIVERGE joins cluster names with '/'."""
    return f"{group_a}/{group_b}"


def normalise(
    results: pd.DataFrame,
    family: str,
    partition: str,
    group_a: str,
    group_b: str,
    test: str,
) -> pd.DataFrame:
    """Map a DIVERGE results frame onto the SITES schema."""
    label = column_label(group_a, group_b)
    if label in results.columns:
        posterior = results[label]
    elif results.shape[1] == 1:
        posterior = results.iloc[:, 0]
    else:
        raise KeyError(
            f"expected a column {label!r} in DIVERGE results, got {list(results.columns)}"
        )

    out = pd.DataFrame(
        {
            "family": family,
            "partition": partition,
            "group_a": group_a,
            "group_b": group_b,
            "column": pd.Series(results.index, dtype="Int64"),
            "test": test,
            "posterior": pd.to_numeric(posterior.to_numpy(), errors="coerce"),
        }
    )
    return SITES.validate(out)


def normalise_summary(
    summary: pd.DataFrame,
    family: str,
    partition: str,
    group_a: str,
    group_b: str,
    test: str,
) -> pd.DataFrame:
    """Per-comparison coefficients from `.summary`."""
    label = column_label(group_a, group_b)
    series = summary[label] if label in summary.columns else summary.iloc[:, 0]
    rows = [
        {
            "family": family,
            "partition": partition,
            "group_a": group_a,
            "group_b": group_b,
            "test": test,
            "parameter": str(name),
            "value": pd.to_numeric(value, errors="coerce"),
        }
        for name, value in series.items()
    ]
    return COMPARISONS.validate(pd.DataFrame(rows))


def run_pair(
    alignment: Path,
    tree_a: Path,
    tree_b: Path,
    family: str,
    partition: str,
    group_a: str,
    group_b: str,
    include_type1: bool = True,
    allow_shim: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Run Type-II and, unless disabled, Type-I for one group pair.

    Returns (sites, comparisons, notes); notes record any shim use.
    """
    diverge = _import_diverge()
    names = [group_a, group_b]
    sites, comparisons, notes = [], [], []

    t2 = diverge.Type2(str(alignment), str(tree_a), str(tree_b), cluster_name=names)
    sites.append(normalise(t2.results, family, partition, group_a, group_b, TYPE2))
    comparisons.append(normalise_summary(t2.summary, family, partition, group_a, group_b, TYPE2))

    if include_type1:
        if needs_shim():
            if not allow_shim:
                raise DivergeUpstreamBug(
                    "Type-I unavailable: " + BROKEN_UPSTREAM["Gu99"] + ". "
                    "Pass allow_shim=True or patch upstream."
                )
            if apply_shim():
                notes.append(SHIM_NOTE)
        gu = diverge.Gu99(str(alignment), str(tree_a), str(tree_b), cluster_name=names)
        sites.append(normalise(gu.results, family, partition, group_a, group_b, TYPE1))
        comparisons.append(
            normalise_summary(gu.summary, family, partition, group_a, group_b, TYPE1)
        )

    return (
        pd.concat(sites, ignore_index=True),
        pd.concat(comparisons, ignore_index=True),
        notes,
    )


def run_partition(
    alignment: Path,
    trees: dict[str, Path],
    partition: Partition,
    family: str,
    skip_underpowered: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    """Run every group pair; skipped pairs are returned, not dropped."""
    skipped: list[str] = []
    notes: list[str] = []
    weak = set(partition.underpowered_groups())
    site_frames, comparison_frames = [], []

    for a, b in partition.pairs():
        if a not in trees or b not in trees:
            skipped.append(f"{a}|{b}: missing conformant tree")
            continue
        if skip_underpowered and (a in weak or b in weak):
            skipped.append(f"{a}|{b}: underpowered")
            continue
        s, c, pair_notes = run_pair(alignment, trees[a], trees[b], family, partition.name, a, b)
        site_frames.append(s)
        comparison_frames.append(c)
        notes.extend(pair_notes)

    empty_sites = SITES.validate(pd.DataFrame(columns=[c.name for c in SITES.columns]))
    empty_comp = COMPARISONS.validate(pd.DataFrame(columns=[c.name for c in COMPARISONS.columns]))
    return (
        pd.concat(site_frames, ignore_index=True) if site_frames else empty_sites,
        pd.concat(comparison_frames, ignore_index=True) if comparison_frames else empty_comp,
        skipped,
        sorted(set(notes)),
    )
