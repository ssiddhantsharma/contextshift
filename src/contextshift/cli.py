"""Command line surface."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from . import __version__, join, report, stats
from .partition import Partition, Provenance, adjusted_rand_index
from .schema import ALL_SCHEMAS
from .stages import (
    conserve,
    derep,
    diverge,
    families,
    loci,
    mapping,
    profiles,
    structures,
    toolchain_report,
)

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


@app.command("loci")
def loci_cmd(
    operons: Path = typer.Argument(..., help="cas_operons.tab from a typer"),
    out: Path = typer.Argument(..., help="hits parquet"),
    partition_out: Path = typer.Option(None, help="also write a subtype partition JSON"),
    source: str = typer.Option("unknown", help="typer name and version"),
    scheme: str = typer.Option("unknown", help="classification scheme and version"),
) -> None:
    """Parse typed loci into per-gene hits, and optionally a subtype partition."""
    hits = loci.to_hits(loci.read_operons(operons))
    hits.to_parquet(out, index=False)
    typer.echo(f"{len(hits)} hits from {hits['locus_id'].nunique()} loci -> {out}")
    if partition_out:
        p, dropped = loci.subtype_partition(hits, source, scheme)
        p.to_json(partition_out)
        typer.echo(f"partition: {p.counts()} -> {partition_out}")
        if dropped:
            typer.echo(f"dropped {len(dropped)} members with an unusable subtype")


@app.command("families")
def families_cmd(
    hits: Path = typer.Argument(..., help="hits parquet from `loci`"),
    out: Path = typer.Argument(..., help="members parquet"),
    fusion_policy: str = typer.Option(families.FLAG, help="flag | exclude"),
) -> None:
    """Group hits into families, handling fused ORFs and paralogous copies."""
    members, rep = families.build(pd.read_parquet(hits), fusion_policy=fusion_policy)
    members.to_parquet(out, index=False)
    typer.echo(rep.as_text())
    typer.echo(f"-> {out}")


@app.command("derep")
def derep_cmd(
    fasta: Path = typer.Argument(..., help="protein FASTA"),
    out: Path = typer.Argument(..., help="derep parquet"),
    identity: float = typer.Option(0.90, help="MMseqs2 --min-seq-id"),
    threads: int = typer.Option(4),
) -> None:
    """Cluster sequences, keeping cluster size as a weight."""
    d = derep.cluster(fasta, identity=identity, threads=threads)
    d.to_parquet(out, index=False)
    n = d["cluster_id"].nunique()
    typer.echo(f"{len(d)} sequences -> {n} clusters ({len(d)/n:.2f}x redundancy) -> {out}")


@app.command("conserve")
def conserve_cmd(
    alignment: Path = typer.Argument(...),
    out: Path = typer.Argument(..., help="conservation parquet"),
    family: str = typer.Option(...),
    scope: str = typer.Option(conserve.SCOPE_ALL, help="'all' or a group name"),
) -> None:
    """Per-column conservation (Jensen-Shannon; lower means more constrained)."""
    d = conserve.jensen_shannon(alignment, family, scope)
    d.to_parquet(out, index=False)
    typer.echo(f"{len(d)} columns, {int(d['rate'].isna().sum())} too gappy to score -> {out}")


@app.command("diverge")
def diverge_cmd(
    alignment: Path = typer.Argument(...),
    tree_a: Path = typer.Argument(...),
    tree_b: Path = typer.Argument(...),
    out: Path = typer.Argument(..., help="sites parquet"),
    family: str = typer.Option(...),
    group_a: str = typer.Option("A"),
    group_b: str = typer.Option("B"),
    comparisons_out: Path = typer.Option(None, help="per-comparison coefficients"),
) -> None:
    """Type-I and Type-II divergence for one group pair."""
    sites, comps, notes = diverge.run_pair(
        alignment, tree_a, tree_b, family, "supplied", group_a, group_b)
    sites.to_parquet(out, index=False)
    typer.echo(f"{len(sites)} site rows over {sites['column'].nunique()} columns -> {out}")
    for n in notes:
        typer.echo(f"note: {n}")
    if comparisons_out:
        comps.to_parquet(comparisons_out, index=False)


@app.command("typers")
def typers_cmd(
    a: Path = typer.Argument(..., help="first typer output"),
    b: Path = typer.Argument(..., help="second typer output"),
    format_a: str = typer.Option("cctyper", help="cctyper | defensefinder | padloc"),
    format_b: str = typer.Option("defensefinder", help="cctyper | defensefinder | padloc"),
    system: str = typer.Option(None, help="regex selecting one system family"),
    coarse: bool = typer.Option(False, help="drop extra subtype granularity before comparing"),
) -> None:
    """Agreement between two typers' partitions, as an adjusted Rand index."""
    readers = {
        "defensefinder": lambda p: loci.parse_defensefinder(p, system=system),
        "padloc": lambda p: loci.parse_padloc(p, system=system),
        "cctyper": lambda p: loci.to_hits(loci.read_operons(p)),
    }
    for name in (format_a, format_b):
        if name not in readers:
            raise typer.BadParameter(f"unknown format {name!r}; choose from {sorted(readers)}")

    frames = []
    for path, fmt in ((a, format_a), (b, format_b)):
        df = readers[fmt](path)
        if coarse:
            df = df.assign(subtype=df["subtype"].map(loci.coarse_subtype))
        frames.append(loci.partition_from(df, "subtype", fmt, "as-reported"))

    pa, pb = frames
    shared = set(pa.labels) & set(pb.labels)
    typer.echo(f"{format_a}: {len(pa.labels)} members, {len(pa.group_names)} groups")
    typer.echo(f"{format_b}: {len(pb.labels)} members, {len(pb.group_names)} groups")
    typer.echo(f"shared members: {len(shared)}")
    typer.echo(f"adjusted rand : {adjusted_rand_index(pa, pb):.4f}")
    disagree = [m for m in sorted(shared) if pa.labels[m] != pb.labels[m]]
    typer.echo(f"disagreeing   : {len(disagree)}")
    for m in disagree[:10]:
        typer.echo(f"  {m}: {format_a}={pa.labels[m]}  {format_b}={pb.labels[m]}")


@app.command("profiles")
def profiles_cmd(
    fasta: Path = typer.Argument(...),
    out: Path = typer.Argument(..., help="assignment parquet"),
    pfam: list[str] = typer.Option(..., "--pfam", help="Pfam accession, repeatable"),
    hmmdir: Path = typer.Option(Path("hmms"), help="where to cache downloaded HMMs"),
    evalue: float = typer.Option(1e-5),
    min_margin: float = typer.Option(10.0, help="bits below which a call is ambiguous"),
) -> None:
    """Assign family membership by profile hit rather than by label."""
    paths = [profiles.fetch_pfam(acc, hmmdir) for acc in pfam]
    combined = profiles.combine(paths, Path(hmmdir) / "combined.hmm")
    hits = profiles.scan(fasta, combined, evalue=evalue)
    assigned = profiles.assign(hits)
    assigned.to_parquet(out, index=False)
    typer.echo(f"{len(assigned)} sequences assigned -> {out}")
    typer.echo(assigned["family"].value_counts().to_string())
    amb = profiles.ambiguous(assigned, min_margin=min_margin)
    if len(amb):
        typer.echo(f"ambiguous (margin < {min_margin} bits): {len(amb)}")


@app.command("structures")
def structures_cmd(
    accessions: Path = typer.Argument(..., help="one UniProt accession per line"),
    outdir: Path = typer.Argument(..., help="where to write mmCIF files"),
    report_to: Path = typer.Option(None, help="write the model table here"),
) -> None:
    """Fetch AlphaFold models, gating on confidence."""
    accs = [a.strip() for a in Path(accessions).read_text().split() if a.strip()]
    df, notes = structures.fetch_many(accs, outdir)
    typer.echo(f"{len(df)} models, {int(df['usable'].sum()) if len(df) else 0} usable")
    for n in notes[:20]:
        typer.echo(f"  {n}")
    if report_to:
        df.to_parquet(report_to, index=False)


@app.command("map")
def map_cmd(
    alignment: Path = typer.Argument(...),
    structure: Path = typer.Argument(..., help="mmCIF"),
    chain: str = typer.Argument(...),
    reference_id: str = typer.Argument(..., help="sequence id in the alignment"),
    out: Path = typer.Argument(..., help="mapping parquet"),
    family: str = typer.Option(...),
) -> None:
    """Map alignment columns onto author residue numbers."""
    residues = mapping.read_chain(structure, chain)
    aligned = mapping.read_alignment(alignment)[reference_id]
    m, mismatches = mapping.map_columns(aligned, residues, family, reference_id)
    m.to_parquet(out, index=False)
    typer.echo(f"{len(m)} columns mapped to {chain}:{min(m['resnum'])}..{max(m['resnum'])} -> {out}")


@app.command("report")
def report_cmd(
    classified: Path = typer.Argument(..., help="classified parquet from `classify`"),
    outdir: Path = typer.Argument(...),
) -> None:
    """Write figures and a text summary."""
    df = pd.read_parquet(classified)
    run = report.RunSummary()
    run.stage("columns", len(df))
    for p in report.write(outdir, run, df):
        typer.echo(f"  {p}")


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
