# contextshift

Find the residues that distinguish one group of a protein family from another,
where the grouping is supplied from outside the phylogeny.

## Why

Functional-divergence methods normally cut a tree into subfamilies and test
those. That works only when the grouping you care about is monophyletic. In
many prokaryotic families it is not: the meaningful grouping comes from a
profile assignment, a gene neighbourhood, or membership of a larger system, and
it cuts across the tree.

`contextshift` takes the grouping as an input and records where it came from,
so a stale or disputed classification is visible rather than silent.

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

## Pipeline

```
typed loci ─families─▶ members ─derep─▶ weights
                                          │
                                       align ─▶ MSA ─tree─▶ per-group Newick
                                                  │               │
                                              conserve        diverge
                                                  └── classify ───┘
                                                        │
                                              map ─▶ sites × class × structure
                                                        │
                                                     report
```

Every stage boundary is a declared table (`contextshift schemas`), so you can
enter anywhere with your own data and be told at once if it does not fit.

## Check power before anything expensive

```bash
contextshift make-partition labels.tsv groups.json \
  --name context --source mytyper-1.2 --scheme "some-classification-2025"
contextshift power groups.json
```

Dereplication keeps cluster size as a weight, and group size is reported as
summed weights, so a hundred copies of one sequence count once. A group that
looks large and is not gets flagged before you spend anything on it.

## Design rules

- The library is target-agnostic. A biological finding in library source fails
  a test; naming a tool it drives does not.
- Reference data from a prior study never replaces a computed result. It
  produces a comparison.
- Multiple-testing correction is global by default, across family × pair × column.
- Skipped comparisons are returned, never silently omitted.

## Verification

- [`METHODS.md`](METHODS.md) — the source for every method, and why each
  substitution was made.
- [`VERIFICATION.md`](VERIFICATION.md) — every assumption about an external
  tool, and whether it was measured, documented, or is still unchecked.

Two properties are asserted against DIVERGE itself: this library runs group
pairs one at a time, and that is only sound because a pairwise estimate is
byte-identical to extracting the same pair from a multi-cluster run
(`tests/test_validation.py`).

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
| [biotite](https://www.biotite-python.org/) | Kunzmann & Hamacher 2018, *BMC Bioinformatics* 19:346 |

Conservation is computed in-library by Jensen-Shannon divergence
(Capra & Singh 2007, *Bioinformatics* 23:1875) rather than by shelling out;
Rate4Site (Pupko et al. 2002) is kept as an alternative.

## Status

Early. Every stage has producing code and every schema a producer. 142 tests,
including live runs of DIVERGE, MAFFT and MMseqs2. The MEME and Rate4Site
adapters are unrun and optional.

Type-I divergence needs a one-line fix for a defect in DIVERGE 4.1.0
([upstream PR](https://github.com/zjupgx/diverge4/pull/8)). The library applies
it as a shim and records that it did, so a Type-I result is never mistaken for
one from stock upstream.

## License

MIT
