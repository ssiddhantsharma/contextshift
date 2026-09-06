"""Typed loci into family hits and a subtype partition.

Three typers, three formats:

CCTyper `cas_operons.tab` holds Python lists written straight into TSV cells by
`to_csv`, so they arrive as `"['a', 'b']"` and are parallel per gene.

DefenseFinder `defense_finder_systems.tsv` names a system's type and subtype
separately, the subtype as e.g. `CAS_Class1-Subtype-I-E`. Its
`protein_in_syst` and `name_of_profiles_in_sys` are each sorted alphabetically
and so do NOT correspond position by position; the gene table is the only safe
source for that mapping.

Both typers detect many system families, so callers pass a `system` pattern to
select the one under study.

PADLOC `_padloc.csv` is one row per protein already, with dots in the column
names.
"""

from __future__ import annotations

import ast
import re
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


# --- DefenseFinder ---------------------------------------------------------

DF_SUBTYPE = re.compile(r"Subtype-([IVX]+-[A-Z]\d?)", re.I)


def parse_defensefinder(systems: Path, system: str | None = None) -> pd.DataFrame:
    """One row per (system, protein) from defense_finder_systems.tsv.

    `protein_in_syst` is expanded, but never zipped against
    `name_of_profiles_in_sys`: the file documents both as independently
    alphabetised, so pairing them by position silently mislabels genes.
    """
    df = pd.read_csv(systems, sep="\t")
    missing = {"sys_id", "type", "subtype", "protein_in_syst"} - set(df.columns)
    if missing:
        raise LocusError(f"{systems}: missing columns {sorted(missing)}")
    if system:
        keep = df["type"].astype(str).str.contains(system, case=False, regex=True) | df[
            "subtype"].astype(str).str.contains(system, case=False, regex=True)
        df = df[keep]

    rows = []
    for r in df.itertuples(index=False):
        m = DF_SUBTYPE.search(str(r.subtype))
        for protein in str(r.protein_in_syst).split(","):
            protein = protein.strip()
            if protein:
                rows.append({
                    "member_id": protein,
                    "locus_id": str(r.sys_id),
                    "subtype": m.group(1).upper() if m else str(r.subtype),
                    "raw_subtype": str(r.subtype),
                })
    return pd.DataFrame(rows, columns=["member_id", "locus_id", "subtype", "raw_subtype"])


# --- PADLOC ----------------------------------------------------------------

# PADLOC names systems `cas_type_I-E`, and subdivides some subtypes that other
# typers do not: I-B1/I-B2, I-F1/I-F2/I-F3. It also has a `cas_adaptation`
# system with no subtype at all. Both are reported, not silently coerced.
PADLOC_SUBTYPE = re.compile(r"type[_\s]*([IVX]+-[A-Z]\d?)", re.I)


def parse_padloc(csv: Path, system: str | None = None) -> pd.DataFrame:
    """One row per protein from a `_padloc.csv`. Column names contain dots."""
    df = pd.read_csv(csv)
    missing = {"system", "target.name", "system.number", "seqid"} - set(df.columns)
    if missing:
        raise LocusError(f"{csv}: missing columns {sorted(missing)}")
    if system:
        df = df[df["system"].astype(str).str.contains(system, case=False, na=False, regex=True)]

    subtype = df["system"].astype(str).map(
        lambda v: (m.group(1).upper() if (m := PADLOC_SUBTYPE.search(v)) else v)
    )
    return pd.DataFrame({
        "member_id": df["target.name"].astype(str),
        "locus_id": df["seqid"].astype(str) + "@" + df["system.number"].astype(str),
        "subtype": subtype,
        "raw_subtype": df["system"].astype(str),
    }).reset_index(drop=True)


def partition_from(labels: pd.DataFrame, name: str, source: str, scheme: str) -> Partition:
    """Build a Partition from any typer's (member_id, subtype) table."""
    df = labels[~labels["subtype"].isin(AMBIGUOUS)].drop_duplicates("member_id")
    if df.empty:
        raise LocusError("no members with a usable subtype")
    return Partition(
        name=name,
        labels=dict(zip(df["member_id"], df["subtype"], strict=True)),
        provenance=Provenance(source=source, scheme=scheme),
    )


def coarse_subtype(subtype: str) -> str:
    """Drop a typer's extra granularity, e.g. PADLOC's I-B1 -> I-B.

    Only for comparing partitions whose schemes differ in depth. It loses
    information, so it is never applied to the labels an analysis runs on.
    """
    m = re.fullmatch(r"([IVX]+-[A-Z])\d", subtype, re.I)
    return m.group(1).upper() if m else subtype
