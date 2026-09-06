"""Group typed gene hits into families, one row per (ORF, gene)."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from ..schema import MEMBERS

FLAG = "flag"
EXCLUDE = "exclude"
POLICIES = (FLAG, EXCLUDE)


@dataclass
class FamilyReport:
    n_hits: int = 0
    n_members: int = 0
    n_fusion_orfs: int = 0
    fusion_genes: dict[str, int] = field(default_factory=dict)
    n_excluded: int = 0
    n_multicopy_loci: int = 0
    max_copies: int = 0
    families: dict[str, int] = field(default_factory=dict)

    def as_text(self) -> str:
        lines = [
            f"hits            {self.n_hits}",
            f"members kept    {self.n_members}",
            f"fusion ORFs     {self.n_fusion_orfs}",
            f"excluded        {self.n_excluded}",
            f"multi-copy loci {self.n_multicopy_loci} (max {self.max_copies} copies)",
        ]
        for combo, n in sorted(self.fusion_genes.items(), key=lambda kv: -kv[1]):
            lines.append(f"  fused: {combo}  x{n}")
        for fam, n in sorted(self.families.items(), key=lambda kv: -kv[1]):
            lines.append(f"  family {fam}: {n}")
        return "\n".join(lines)


def build(
    hits: pd.DataFrame,
    fusion_policy: str = FLAG,
    order_by: str = "start",
) -> tuple[pd.DataFrame, FamilyReport]:
    """Turn gene hits into a MEMBERS table.

    Needs member_id, genome_id, gene, locus_id; uses start/end/strand and
    hit_evalue when present. Copies of one gene within a locus are numbered
    from 1 in `order_by` order.
    """
    if fusion_policy not in POLICIES:
        raise ValueError(f"fusion_policy must be one of {POLICIES}")

    required = {"member_id", "genome_id", "gene", "locus_id"}
    missing = required - set(hits.columns)
    if missing:
        raise ValueError(f"hits missing columns: {sorted(missing)}")

    df = hits.copy()
    report = FamilyReport(n_hits=len(df))

    genes_per_orf = df.groupby("member_id")["gene"].transform("nunique")
    df["is_fusion"] = genes_per_orf > 1

    fused = df[df["is_fusion"]]
    report.n_fusion_orfs = int(fused["member_id"].nunique())
    for _, grp in fused.groupby("member_id"):
        combo = "+".join(sorted(set(grp["gene"])))
        report.fusion_genes[combo] = report.fusion_genes.get(combo, 0) + 1

    if fusion_policy == EXCLUDE:
        report.n_excluded = int(df["is_fusion"].sum())
        df = df[~df["is_fusion"]].copy()

    sort_cols = [c for c in (order_by, "hit_evalue") if c in df.columns]
    if sort_cols:
        df = df.sort_values(["locus_id", "gene", *sort_cols], kind="stable")
    df["copy_index"] = df.groupby(["locus_id", "gene"]).cumcount() + 1

    copies = df.groupby(["locus_id", "gene"])["copy_index"].max()
    report.n_multicopy_loci = int((copies > 1).sum())
    report.max_copies = int(copies.max()) if len(copies) else 0

    df["family"] = df["gene"]
    report.families = df["family"].value_counts().to_dict()
    report.n_members = len(df)

    keep = [c for c in MEMBERS.required_names + [
        "contig", "start", "end", "strand", "is_fusion", "copy_index", "hit_evalue"
    ] if c in df.columns]
    out = df[keep].reset_index(drop=True)

    # a fused ORF appears once per gene identity
    out = out.drop_duplicates(subset=["member_id", "family"])
    return MEMBERS.validate(out.assign(member_id=out["member_id"].astype("string"))), report


def family_key(row) -> str:
    """Family name, suffixed when it is not the first copy in a locus."""
    idx = int(row.copy_index) if pd.notna(row.copy_index) else 1
    return f"{row.family}-{idx}" if idx > 1 else str(row.family)


def split_by_partition(members: pd.DataFrame, labels: dict[str, str]) -> dict[str, pd.DataFrame]:
    """Split members by label, keeping unlabelled members visible."""
    df = members.copy()
    df["_group"] = df["member_id"].map(labels).fillna("unlabelled")
    return {g: sub.drop(columns="_group") for g, sub in df.groupby("_group")}
