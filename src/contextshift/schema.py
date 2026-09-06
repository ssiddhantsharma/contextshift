"""Stage boundary contracts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


class SchemaError(ValueError):
    pass


@dataclass(frozen=True)
class Column:
    name: str
    dtype: str
    required: bool = True
    nullable: bool = False


@dataclass(frozen=True)
class TableSchema:
    name: str
    columns: tuple[Column, ...]
    key: tuple[str, ...] = ()

    @property
    def required_names(self) -> list[str]:
        return [c.name for c in self.columns if c.required]

    def validate(self, df: pd.DataFrame) -> pd.DataFrame:
        missing = [c for c in self.required_names if c not in df.columns]
        if missing:
            raise SchemaError(f"{self.name}: missing required columns {missing}")

        for col in self.columns:
            if col.name not in df.columns:
                continue
            if not col.nullable and df[col.name].isna().any():
                n = int(df[col.name].isna().sum())
                raise SchemaError(f"{self.name}: column {col.name!r} has {n} null values")
            try:
                df[col.name] = df[col.name].astype(col.dtype)
            except (TypeError, ValueError) as exc:
                raise SchemaError(
                    f"{self.name}: column {col.name!r} not castable to {col.dtype}: {exc}"
                ) from exc

        if self.key:
            dup = df.duplicated(subset=list(self.key), keep=False)
            if dup.any():
                raise SchemaError(
                    f"{self.name}: {int(dup.sum())} rows duplicate key {self.key}"
                )
        return df

    def read(self, path: Path) -> pd.DataFrame:
        return self.validate(pd.read_parquet(path))

    def write(self, df: pd.DataFrame, path: Path) -> Path:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        self.validate(df).to_parquet(path, index=False)
        return Path(path)


MEMBERS = TableSchema(
    name="members",
    columns=(
        Column("member_id", "string"),
        Column("genome_id", "string"),
        Column("family", "string"),
        Column("contig", "string", required=False, nullable=True),
        Column("start", "Int64", required=False, nullable=True),
        Column("end", "Int64", required=False, nullable=True),
        Column("strand", "string", required=False, nullable=True),
        Column("is_fusion", "boolean", required=False, nullable=True),
        Column("copy_index", "Int64", required=False, nullable=True),
        Column("hit_evalue", "float64", required=False, nullable=True),
    ),
    # a fused ORF is one member_id under two families
    key=("member_id", "family"),
)

DEREP = TableSchema(
    name="derep",
    columns=(
        Column("member_id", "string"),
        Column("representative_id", "string"),
        Column("cluster_id", "string"),
        Column("cluster_size", "Int64"),
        Column("weight", "float64"),
        Column("level", "string"),
    ),
    key=("member_id", "level"),
)

SITES = TableSchema(
    name="sites",
    columns=(
        Column("family", "string"),
        Column("partition", "string"),
        Column("group_a", "string"),
        Column("group_b", "string"),
        # 0-based; DIVERGE returns a sparse set of kept positions
        Column("column", "Int64"),
        Column("test", "string"),
        # Qk in [0, 1]; not a p-value, never feed to a p-value correction
        Column("posterior", "float64", nullable=True),
        # only for methods that produce p-values
        Column("pvalue", "float64", required=False, nullable=True),
        Column("qvalue", "float64", required=False, nullable=True),
    ),
    key=("family", "partition", "group_a", "group_b", "column", "test"),
)

COMPARISONS = TableSchema(
    name="comparisons",
    columns=(
        Column("family", "string"),
        Column("partition", "string"),
        Column("group_a", "string"),
        Column("group_b", "string"),
        Column("test", "string"),
        # theta, alpha, SE: properties of the group pair, not of a site
        Column("parameter", "string"),
        Column("value", "float64", nullable=True),
    ),
    key=("family", "partition", "group_a", "group_b", "test", "parameter"),
)

CONSERVATION = TableSchema(
    name="conservation",
    columns=(
        Column("family", "string"),
        Column("column", "Int64"),
        Column("scope", "string"),
        Column("rate", "float64", nullable=True),
        Column("score", "float64", required=False, nullable=True),
        Column("occupancy", "float64", required=False, nullable=True),
    ),
    key=("family", "column", "scope"),
)

NEIGHBOURS = TableSchema(
    name="neighbours",
    columns=(
        Column("member_id", "string"),
        Column("offset", "Int64"),
        Column("neighbour_id", "string"),
        Column("cluster_id", "string", nullable=True),
        Column("annotation", "string", required=False, nullable=True),
        Column("conservation_rank", "Int64", required=False, nullable=True),
    ),
    key=("member_id", "offset"),
)

MOTIFS = TableSchema(
    name="motifs",
    columns=(
        Column("family", "string"),
        Column("partition", "string"),
        Column("group", "string"),
        Column("motif_id", "string"),
        Column("consensus", "string"),
        Column("width", "Int64"),
        Column("evalue", "float64", nullable=True),
        Column("n_sites", "Int64"),
    ),
    key=("family", "partition", "group", "motif_id"),
)

MAPPING = TableSchema(
    name="mapping",
    columns=(
        Column("family", "string"),
        Column("column", "Int64"),
        Column("reference_id", "string"),
        Column("resnum", "Int64", nullable=True),
        Column("resname", "string", required=False, nullable=True),
    ),
    key=("family", "column", "reference_id"),
)

ALL_SCHEMAS: dict[str, TableSchema] = {
    s.name: s
    for s in (MEMBERS, DEREP, SITES, COMPARISONS, CONSERVATION, NEIGHBOURS, MOTIFS, MAPPING)
}
