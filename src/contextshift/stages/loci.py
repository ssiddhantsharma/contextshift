"""Typed loci into family hits and a subtype partition.

Parses CCTyper `cas_operons.tab`. Its `Genes`, `Positions`, `E-values`,
`CoverageSeq` and `CoverageHMM` columns hold Python lists written straight into
TSV cells by `to_csv`, so they arrive as `"['a', 'b']"` and are parallel: one
entry per gene in the locus.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pandas as pd

from ..partition import Partition, Provenance

LIST_COLUMNS = ("Genes", "Positions", "E-values", "CoverageSeq", "CoverageHMM")
REQUIRED = ("Contig", "Operon", "Prediction", "Genes", "Positions")
AMBIGUOUS = ("False", "Ambiguous", "Unknown")


class LocusError(ValueError):
    pass


def _as_list(value):
    if isinstance(value, list):
        return value
    if pd.isna(value):
        return []
    parsed = ast.literal_eval(str(value))
    return list(parsed) if isinstance(parsed, (list, tuple)) else [parsed]


def read_operons(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise LocusError(f"{path}: missing columns {missing}")
    for col in LIST_COLUMNS:
        if col in df.columns:
            df[col] = df[col].map(_as_list)
    return df


def to_hits(operons: pd.DataFrame, genome_id: str | None = None) -> pd.DataFrame:
    """Explode parallel per-gene lists into one row per gene."""
    rows = []
    for r in operons.itertuples(index=False):
        genes, positions = list(r.Genes), list(r.Positions)
        if len(genes) != len(positions):
            raise LocusError(
                f"operon {r.Operon}: {len(genes)} genes but {len(positions)} positions"
            )
        for gene, pos in zip(genes, positions, strict=True):
            rows.append(
                {
                    "member_id": str(pos),
                    "genome_id": genome_id or str(r.Contig),
                    "contig": str(r.Contig),
                    "locus_id": str(r.Operon),
                    "gene": str(gene),
                    "subtype": str(r.Prediction),
                    "start": int(r.Start) if hasattr(r, "Start") else None,
                }
            )
    columns = ["member_id", "genome_id", "contig", "locus_id", "gene", "subtype", "start"]
    return pd.DataFrame(rows, columns=columns)


def subtype_partition(
    hits: pd.DataFrame,
    typer_version: str,
    scheme: str,
    drop_ambiguous: bool = True,
) -> tuple[Partition, list[str]]:
    """Partition members by locus subtype.

    Returns the partition and the members dropped for an unusable label, so a
    dropped locus is counted rather than disappearing.
    """
    df = hits.copy()
    dropped = sorted(df[df["subtype"].isin(AMBIGUOUS)]["member_id"]) if drop_ambiguous else []
    if drop_ambiguous:
        df = df[~df["subtype"].isin(AMBIGUOUS)]
    if df.empty:
        raise LocusError("no members left after dropping ambiguous subtypes")

    return (
        Partition(
            name="subtype",
            labels=dict(zip(df["member_id"], df["subtype"], strict=True)),
            provenance=Provenance(source=typer_version, scheme=scheme),
        ),
        dropped,
    )


def gene_presence(operons: pd.DataFrame, gene: str) -> pd.DataFrame:
    """Per subtype, in how many loci a gene occurs. Measures presence, not assumes it."""
    rows = []
    for r in operons.itertuples(index=False):
        rows.append({"subtype": str(r.Prediction), "present": gene.lower() in
                     [str(g).lower() for g in r.Genes]})
    df = pd.DataFrame(rows)
    out = df.groupby("subtype")["present"].agg(["sum", "count"]).reset_index()
    out.columns = ["subtype", "n_loci_with_gene", "n_loci"]
    out["fraction"] = out["n_loci_with_gene"] / out["n_loci"]
    return out.sort_values("fraction", ascending=False).reset_index(drop=True)
