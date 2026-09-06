"""Multiple sequence alignment.

Structure-guided alignment is the default because site-level tests inherit
every alignment error, and these families are divergent enough that sequence
alignment alone puts non-homologous columns together.
"""

from __future__ import annotations

from pathlib import Path

from ._external import FOLDMASON, MAFFT, run

SEQUENCE = "sequence"
STRUCTURE = "structure"


def mafft(fasta: Path, out: Path, threads: int = 4, accurate: bool = True) -> Path:
    MAFFT.require()
    args = ["mafft", "--thread", str(threads), "--anysymbol"]
    args += ["--maxiterate", "1000", "--localpair"] if accurate else ["--auto"]
    args.append(str(fasta))
    result = run(args)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(result.stdout)
    return Path(out)


def foldmason(structure_dir: Path, out_prefix: Path, threads: int = 4) -> Path:
    """Structure-guided MSA from predicted or experimental structures."""
    FOLDMASON.require()
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            "foldmason", "easy-msa", str(structure_dir), str(out_prefix),
            str(Path(out_prefix).parent / "foldmason_tmp"),
            "--threads", str(threads),
        ]
    )
    return Path(str(out_prefix) + "_aa.fa")


def occupancy(alignment: Path) -> list[float]:
    """Per-column fraction of non-gap residues."""
    from Bio import AlignIO

    aln = AlignIO.read(str(alignment), "fasta")
    n = len(aln)
    return [
        sum(1 for rec in aln if rec.seq[i] not in "-.") / n
        for i in range(aln.get_alignment_length())
    ]
