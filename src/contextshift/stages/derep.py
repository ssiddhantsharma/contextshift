"""Dereplication with cluster size retained as a weight.

Weights let effective sample size stay computable downstream.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd

from ..schema import DEREP
from ._external import MMSEQS, run

GENOME = "genome"
SEQUENCE = "sequence"


def cluster(
    fasta: Path,
    identity: float,
    coverage: float = 0.8,
    level: str = SEQUENCE,
    threads: int = 4,
    workdir: Path | None = None,
) -> pd.DataFrame:
    """Cluster with MMseqs2 and return a DEREP table."""
    MMSEQS.require()
    tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="contextshift-derep-"))
    tmp.mkdir(parents=True, exist_ok=True)
    prefix = tmp / "clu"

    run(
        [
            "mmseqs", "easy-cluster", str(fasta), str(prefix), str(tmp / "mmseqs_tmp"),
            "--min-seq-id", str(identity),
            "-c", str(coverage),
            "--cov-mode", "0",
            "--threads", str(threads),
        ]
    )
    return read_mmseqs_tsv(prefix.with_name(prefix.name + "_cluster.tsv"), level=level)


def read_mmseqs_tsv(path: Path, level: str = SEQUENCE) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", header=None, names=["representative_id", "member_id"])
    return build(raw, level=level)


def build(pairs: pd.DataFrame, level: str = SEQUENCE) -> pd.DataFrame:
    """Turn (representative_id, member_id) pairs into a weighted DEREP table."""
    df = pairs.copy()
    sizes = df.groupby("representative_id")["member_id"].transform("size")
    df["cluster_id"] = df["representative_id"]
    df["cluster_size"] = sizes.astype("Int64")
    df["weight"] = 1.0 / sizes
    df["level"] = level
    return DEREP.validate(
        df[["member_id", "representative_id", "cluster_id", "cluster_size", "weight", "level"]]
    )


def representative_weights(derep: pd.DataFrame) -> dict[str, float]:
    """One unit of independent evidence per cluster."""
    df = DEREP.validate(derep.copy())
    reps = df.drop_duplicates(subset=["representative_id"])
    return {str(r.representative_id): 1.0 for r in reps.itertuples(index=False)}


def member_weights(derep: pd.DataFrame) -> dict[str, float]:
    """Cluster mass shared across its members."""
    df = DEREP.validate(derep.copy())
    return {str(r.member_id): float(r.weight) for r in df.itertuples(index=False)}
