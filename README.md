# contextshift

[![ci](https://github.com/ssiddhantsharma/contextshift/actions/workflows/ci.yml/badge.svg)](https://github.com/ssiddhantsharma/contextshift/actions/workflows/ci.yml)
[![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](pyproject.toml)

Find the residues that distinguish one group of a protein family from another,
where the grouping is supplied from outside the phylogeny.

## Why

Functional-divergence methods cut a tree into subfamilies and test those. That
works only when the grouping you care about is monophyletic. Often it is not:
the meaningful grouping comes from a profile assignment, a gene neighbourhood,
or membership of a larger system, and it cuts across the tree.

`contextshift` takes the grouping as input and records where it came from, so a
stale or disputed classification stays visible.

## The result is a 2×2

Divergence alone cannot say whether a site matters. Conservation alone cannot
see a between-group shift, because it averages over the groups. The cross of
the two is the result.

|  | not divergent | divergent |
|---|---|---|
| conserved everywhere | `core` — catalytic; the negative control | `suspect` — usually misalignment |
| conserved within groups, differs between | — | **`determinant`** |
| constrained in one group only | — | `relaxed` |
| variable everywhere | `variable` | `variable` |

A column whose inputs are missing is emitted with `status` set. Nothing is
dropped quietly; low-occupancy columns are flagged, not filtered.

## Install

```bash
uv venv --python 3.12
uv pip install -e ".[dev]"
contextshift doctor
```

`doctor` reports which optional tools are on PATH. All are checked at call
time, so the library imports without them.

DIVERGE is not a declared dependency: it cannot be installed from PyPI, because
its sdist reads a `requirements.txt` it does not ship and its default source
tree is Windows-only. Build it from a clone using `src_linux`; see
[VERIFICATION.md](VERIFICATION.md).

The `Dockerfile` does all of that for you, and brings MMseqs2, MAFFT, IQ-TREE,
MEME, pyhmmer and DefenseFinder with it. PADLOC is left out of the default
image because its R stack roughly doubles the download; one commented line adds
it:

```bash
docker build -t contextshift .
docker run --rm -v "$PWD:/work" contextshift doctor
```

## Pipeline

```
typed loci ─┬─ families ─▶ members ─derep─▶ weights
            │                                  │
            └─ typers ─▶ agreement          profiles ─▶ family by HMM hit
                         (ARI)                 │
                                            align ─▶ MSA ─tree─▶ per-group Newick
                                                       │               │
                                                   conserve        diverge
                                                       └── classify ───┘
                                                             │
                                    structures ─▶ map ─▶ sites × class × structure
                                                             │
                                                          report
```

Every stage boundary is a declared table (`contextshift schemas`), so you can
enter anywhere with your own data.

```bash
contextshift loci      cas_operons.tab hits.parquet \
                       --partition-out subtype.json --source mytyper-1.2 --scheme scheme-2020
contextshift families  hits.parquet members.parquet
contextshift derep     family.faa derep.parquet --identity 0.90
contextshift conserve  family.aln cons.parquet --family Cas4 --scope all
contextshift diverge   family.aln groupA.nwk groupB.nwk sites.parquet --family Cas4
contextshift fdr       sites.parquet sites_q.parquet
contextshift classify  sites_q.parquet cons.parquet classified.parquet
contextshift map       family.aln ref.cif A ref_id mapping.parquet --family Cas4
contextshift report    classified.parquet out/
```

Worth running before the analysis, not after — whether two typers agree, and
whether family membership is supported by a profile rather than a label:

```bash
contextshift typers   cctyper_out/cas_operons.tab defense_finder_systems.tsv \
                      --format-a cctyper --format-b defensefinder --system CAS
contextshift profiles family.faa assigned.parquet --pfam PF01930 --pfam PF06023
```

## Check power before anything expensive

```bash
contextshift make-partition labels.tsv groups.json \
  --name context --source mytyper-1.2 --scheme "some-classification-2025"
contextshift power groups.json
```

Dereplication keeps cluster size as a weight, and group size is reported as
summed weights, so a hundred copies of one sequence count once. A group that
looks large and is not gets flagged before you spend anything.

## Design rules

- The library is target-agnostic. A biological finding in library source fails
  a test; naming a tool it drives does not. That guard has rejected real code.
- A label is not evidence. Family membership can be assigned by profile hit,
  and a partition can be scored against a second, independent typer.
- Reference data from a prior study never replaces a computed result. It
  produces a comparison.
- Multiple-testing correction is global by default, across family × pair × column.
- Skipped comparisons are returned, never silently omitted.

## Verification

Two stages are checked against published ground truth, using the data
supporting Capra & Singh 2007:

| | |
|---|---|
| catalytic sites score as more constrained | **57 of 57** proteins, mean AUC **0.924** |
| residue-to-ligand distances match published values | **7 of 8** structures exact |

```bash
python benchmarks/capra_singh.py --workdir /tmp/cs   # reproduces both
```

Fixtures ship, so the tests run offline. The eighth structure's difference, and
every assumption this library makes about an external tool, are in
[VERIFICATION.md](VERIFICATION.md). Method sources are in [METHODS.md](METHODS.md).

## Tools

| | |
|---|---|
| [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper) | Russell et al. 2020, *CRISPR J* 3:462 |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | Steinegger & Söding 2017, *Nat Biotechnol* 35:1026 |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | Katoh & Standley 2013, *MBE* 30:772 |
| [FoldMason](https://github.com/steineggerlab/foldmason) | Gilchrist et al. 2024 |
| [IQ-TREE 2](http://www.iqtree.org/) | Minh et al. 2020, *MBE* 37:1530 |
| [DIVERGE v4](https://github.com/zjupgx/diverge4) | Gu 1999, *MBE* 16:1664; Cheng et al. 2025, *MBE* 42:msaf277 |
| [MEME](https://meme-suite.org/) | Bailey & Elkan 1994, *ISMB* 2:28 |
| [FlaGs](https://github.com/GCA-VH-lab/FlaGs) | Saha et al. 2021, *Bioinformatics* 37:1312 |
| [DefenseFinder](https://github.com/mdmparis/defense-finder) | Tesson et al. 2022, *Nat Commun* 13:2561 |
| [PADLOC](https://github.com/padlocbio/padloc) | Payne et al. 2022, *NAR* 50:W541 |
| [pyhmmer](https://github.com/althonos/pyhmmer) | Larralde & Zeller 2023, *Bioinformatics* 39:btad214 |
| [AlphaFold DB](https://alphafold.ebi.ac.uk/) | Varadi et al. 2024, *NAR* 52:D368 |
| [Pfam](https://www.ebi.ac.uk/interpro/) via InterPro | Mistry et al. 2021, *NAR* 49:D412 |
| [biotite](https://www.biotite-python.org/) | Kunzmann & Hamacher 2018, *BMC Bioinformatics* 19:346 |

Conservation is computed in-library by Jensen-Shannon divergence
(Capra & Singh 2007, *Bioinformatics* 23:1875) rather than by shelling out;
Rate4Site (Pupko et al. 2002) is kept as an alternative.

## Status
Type-I divergence needs a one-line fix for a defect in DIVERGE 4.1.0
([upstream PR](https://github.com/zjupgx/diverge4/pull/8)). The library applies
it as a shim and records that it did, so a Type-I result is never mistaken for
one from stock upstream.

## License

MIT
