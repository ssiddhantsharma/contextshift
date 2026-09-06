"""Maximum-likelihood phylogeny with UFBoot supports (>=1000 replicates)."""

from __future__ import annotations

from pathlib import Path

from ._external import IQTREE, run


def iqtree(
    alignment: Path,
    prefix: Path,
    model: str = "MFP",
    ufboot: int = 1000,
    threads: str = "AUTO",
) -> Path:
    IQTREE.require()
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            "iqtree2", "-s", str(alignment), "--prefix", str(prefix),
            "-m", model, "-B", str(ufboot), "-T", str(threads), "--quiet", "-redo",
        ]
    )
    return Path(str(prefix) + ".treefile")


def read_newick(path: Path) -> str:
    return Path(path).read_text().strip()
