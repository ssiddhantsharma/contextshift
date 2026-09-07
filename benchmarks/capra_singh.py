"""Reproduce the published-benchmark numbers quoted in VERIFICATION.md.

Downloads the data supporting Capra & Singh 2007 (Bioinformatics 23:1875) and
runs two checks against it:

  catalytic  does per-site conservation rank Catalytic Site Atlas positions as
             more constrained?
  ligand     do residue-to-ligand distances reproduce the published values?

    python benchmarks/capra_singh.py --workdir /tmp/cs --proteins 60

The repository ships a small fixture from each so the tests run offline; this
script reproduces the full numbers.
"""

from __future__ import annotations

import argparse
import random
import re
import statistics
import tarfile
import urllib.request
from pathlib import Path

DATA_URL = "https://compbio.cs.princeton.edu/conservation/conservation_data.tar.gz"
RCSB = "https://files.rcsb.org/download/{pdb}.cif"


def fetch_data(workdir: Path) -> Path:
    root = workdir / "conservation_data"
    if root.exists():
        return root
    workdir.mkdir(parents=True, exist_ok=True)
    archive = workdir / "conservation_data.tar.gz"
    if not archive.exists():
        print(f"downloading {DATA_URL} (~60 MB)")
        with urllib.request.urlopen(DATA_URL, timeout=600) as r:
            archive.write_bytes(r.read())
    with tarfile.open(archive) as t:
        t.extractall(workdir)
    return root


def catalytic(root: Path, workdir: Path, n_proteins: int, seed: int) -> None:
    from Bio import AlignIO

    from contextshift.stages import conserve

    aln_dir = root / "conservation_alignments" / "csa_hssp"
    pairs = [
        (aln_dir / (c.stem + ".aln"), c)
        for c in sorted((root / "cat_sites").glob("*.cat_sites"))
        if (aln_dir / (c.stem + ".aln")).exists()
    ]
    print(f"\ncatalytic: {len(pairs)} alignment/label pairs available")

    random.seed(seed)
    fa_dir = workdir / "fasta"
    fa_dir.mkdir(exist_ok=True)
    deltas, aucs, cats, nons = [], [], [], []

    for aln_path, cat_path in random.sample(pairs, min(n_proteins, len(pairs))):
        try:
            aln = AlignIO.read(str(aln_path), "clustal")
        except Exception:
            continue
        fa = fa_dir / (aln_path.stem + ".fa")
        AlignIO.write(aln, str(fa), "fasta")
        labels = {
            int(p[0]): int(p[1])
            for p in (line.split() for line in cat_path.read_text().splitlines())
            if len(p) == 2
        }
        d = conserve.jensen_shannon(fa, aln_path.stem, "all").dropna(subset=["rate"])
        d["cat"] = d["column"].map(labels)
        cat, non = d[d["cat"] == 1]["rate"], d[d["cat"] == -1]["rate"]
        if len(cat) < 1 or len(non) < 10:
            continue
        deltas.append(non.mean() - cat.mean())
        aucs.append(sum((non > r).sum() for r in cat) / (len(cat) * len(non)))
        cats.extend(cat)
        nons.extend(non)

    n = len(deltas)
    pooled = sum(sum(1 for x in nons if x > r) for r in cats) / (len(cats) * len(nons))
    print(f"  proteins scored                    {n}")
    print(f"  catalytic more constrained         {sum(d > 0 for d in deltas)}/{n}")
    print(f"  mean rate difference               {statistics.mean(deltas):+.4f}")
    print(f"  mean per-protein AUC               {statistics.mean(aucs):.3f}")
    print(f"  pooled AUC                         {pooled:.3f}")


def ligand(root: Path, workdir: Path, n_structures: int) -> None:
    from contextshift.stages import mapping

    cif_dir = workdir / "cif"
    cif_dir.mkdir(exist_ok=True)
    entries = []
    for f in sorted((root / "lig_distance").rglob("*.dist_to_lig")):
        header = f.read_text(errors="replace").split("\n", 1)[0]
        m = re.match(r"#(\w{4})\s+(\S)\s+(\S+)", header)
        if m and m.group(2) != "_":
            entries.append((f, m.group(1).upper(), m.group(2), m.group(3)))
    print(f"\nligand: {len(entries)} entries with a parseable header")

    exact, checked = 0, 0
    for dist_path, pdb, chain, lig in entries[:n_structures]:
        cif = cif_dir / f"{pdb}.cif"
        if not cif.exists():
            try:
                with urllib.request.urlopen(RCSB.format(pdb=pdb), timeout=120) as r:
                    cif.write_bytes(r.read())
            except Exception as exc:
                print(f"  {pdb}: download failed ({exc})")
                continue
        published = {
            int(p[0]): float(p[2])
            for p in (line.split() for line in dist_path.read_text(errors="replace").splitlines()[1:])
            if len(p) == 3
        }
        try:
            residues = mapping.read_chain(cif, chain)
            mine = mapping.distance_to(cif, chain, residues, target_res_names=(lig,))
        except Exception as exc:
            print(f"  {pdb} {chain} {lig}: {type(exc).__name__}: {str(exc)[:50]}")
            continue
        diffs = [
            abs(published[int(r.resnum)] - r.min_distance)
            for r in mine.itertuples(index=False)
            if int(r.resnum) in published and r.min_distance == r.min_distance
        ]
        if len(diffs) < 20:
            continue
        checked += 1
        worst = max(diffs)
        if worst < 0.01:
            exact += 1
        # deposited hydrogens explain a residual of roughly one bond length
        import biotite.structure as struc  # noqa: F401
        from biotite.structure.io.pdbx import CIFFile, get_structure

        s = get_structure(CIFFile.read(str(cif)), model=1)
        n_h = int((s.element == "H").sum())
        print(f"  {pdb} {chain} {lig:<4} n={len(diffs):>4} max|diff|={worst:.2f} A"
              f"  hydrogens={n_h}")
    print(f"  structures compared                {checked}")
    print(f"  reproduced exactly (<0.01 A)       {exact}/{checked}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", type=Path, default=Path("benchmark_work"))
    ap.add_argument("--proteins", type=int, default=60, help="catalytic proteins to sample")
    ap.add_argument("--structures", type=int, default=8, help="ligand structures to compare")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--only", choices=["catalytic", "ligand"])
    args = ap.parse_args()

    root = fetch_data(args.workdir)
    if args.only in (None, "catalytic"):
        catalytic(root, args.workdir, args.proteins, args.seed)
    if args.only in (None, "ligand"):
        ligand(root, args.workdir, args.structures)


if __name__ == "__main__":
    main()
