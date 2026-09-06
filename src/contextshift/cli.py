"""Command line surface."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from . import __version__, join, stats
from .partition import Partition, Provenance, adjusted_rand_index
from .schema import ALL_SCHEMAS
from .stages import diverge, toolchain_report

app = typer.Typer(add_completion=False, help="Partition-aware functional divergence.")


@app.command()
def version() -> None:
    typer.echo(f"contextshift {__version__}")


@app.command()
def doctor() -> None:
    """Report which external tools are on PATH and whether diverge4 imports."""
    typer.echo(f"contextshift {__version__}")
    typer.echo("external tools:")
    typer.echo(toolchain_report())
    mark = "x" if diverge.available() else " "
    typer.echo(f"  [{mark}] diverge4    (python package)")


@app.command()
def schemas(name: str = typer.Argument(None, help="Schema to describe; omit to list all")) -> None:
    """Print stage table contracts."""
    if name is None:
        for key, schema in ALL_SCHEMAS.items():
            typer.echo(f"{key:14s} key={schema.key} columns={len(schema.columns)}")
        return
    if name not in ALL_SCHEMAS:
        raise typer.BadParameter(f"unknown schema {name!r}; choose from {sorted(ALL_SCHEMAS)}")
    schema = ALL_SCHEMAS[name]
    typer.echo(f"{schema.name}  key={schema.key}")
    for col in schema.columns:
        flags = "required" if col.required else "optional"
        if col.nullable:
            flags += ", nullable"
        typer.echo(f"  {col.name:20s} {col.dtype:10s} ({flags})")


@app.command()
def power(
    partition: Path = typer.Argument(..., help="Partition JSON"),
    min_effective: float = typer.Option(30.0, help="Effective-N floor for a usable group"),
) -> None:
    """Report per-group size and effective size. Run this before anything expensive."""
    p = Partition.from_json(partition)
    typer.echo(f"partition: {p.name}  ({p.provenance.source}, scheme {p.provenance.scheme})")
    typer.echo(f"{'group':30s} {'n':>8s} {'effective':>10s}  usable")
    for row in p.power():
        usable = "yes" if row.effective_n >= min_effective else "NO"
        typer.echo(f"{row.label:30s} {row.n:8d} {row.effective_n:10.1f}  {usable}")
    weak = [r.label for r in p.power() if r.effective_n < min_effective]
    if weak:
        typer.echo(f"\nunderpowered: {', '.join(weak)}")
    for group, caveat in p.qualified_groups().items():
        typer.echo(f"caveat [{group}]: {caveat}")


@app.command()
def concord(
    a: Path = typer.Argument(..., help="Partition JSON"),
    b: Path = typer.Argument(..., help="Partition JSON to compare against"),
) -> None:
    """Adjusted Rand index between two partitions over their shared members."""
    pa, pb = Partition.from_json(a), Partition.from_json(b)
    shared = set(pa.labels) & set(pb.labels)
    typer.echo(f"shared members : {len(shared)}")
    typer.echo(f"{pa.name} groups: {len(pa.group_names)}")
    typer.echo(f"{pb.name} groups: {len(pb.group_names)}")
    typer.echo(f"adjusted rand  : {adjusted_rand_index(pa, pb):.4f}")


@app.command("fdr")
def fdr_cmd(
    sites: Path = typer.Argument(..., help="sites parquet"),
    out: Path = typer.Argument(..., help="output parquet"),
    scope: str = typer.Option(stats.GLOBAL, help="global | family | pair"),
) -> None:
    """Add q-values, corrected globally by default."""
    df = stats.add_qvalues(pd.read_parquet(sites), scope=scope)
    df.to_parquet(out, index=False)
    burden = stats.testing_burden(df)
    typer.echo(f"corrected {len(df)} tests ({scope}) -> {out}")
    typer.echo(burden.to_string(index=False))


@app.command("classify")
def classify_cmd(
    sites: Path = typer.Argument(..., help="sites parquet"),
    conservation: Path = typer.Argument(..., help="conservation parquet"),
    out: Path = typer.Argument(..., help="output parquet"),
    alpha: float = typer.Option(0.05),
    conserved_rate: float = typer.Option(0.5),
    variable_rate: float = typer.Option(1.5),
) -> None:
    """Cross divergence with conservation and classify every column."""
    result = join.classify(
        pd.read_parquet(sites),
        pd.read_parquet(conservation),
        join.JoinThresholds(
            conserved_rate=conserved_rate, variable_rate=variable_rate, alpha=alpha
        ),
    )
    result.to_parquet(out, index=False)
    typer.echo(f"classified {len(result)} columns -> {out}")
    typer.echo(join.summary(result).to_string(index=False))


@app.command("make-partition")
def make_partition(
    tsv: Path = typer.Argument(..., help="member_id<TAB>label[<TAB>weight]"),
    out: Path = typer.Argument(..., help="output JSON"),
    name: str = typer.Option(..., help="Partition name, e.g. context or subtype"),
    source: str = typer.Option(..., help="Where labels came from, e.g. mytyper-1.8.0"),
    scheme: str = typer.Option("unknown", help="Classification scheme and version"),
) -> None:
    """Build a partition from a label table, recording where the labels came from."""
    p = Partition.from_tsv(tsv, name=name, provenance=Provenance(source=source, scheme=scheme))
    p.to_json(out)
    typer.echo(json.dumps({"name": p.name, "members": len(p.labels), "groups": p.counts()}, indent=2))


if __name__ == "__main__":
    app()
