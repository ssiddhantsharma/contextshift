# contextshift

Partition-aware functional divergence for prokaryotic protein families.

Give it a family of proteins and a labelling of that family **from outside the tree** —
genomic context, system subtype, host range — and it reports the sites that discriminate
the labels, with conservation and structural context attached.

## Why the partition comes from outside

Standard functional-divergence workflows cut a phylogeny into subfamilies and test those.
That only works when the grouping you care about is monophyletic. For many prokaryotic
families it is not: the biologically meaningful grouping is assigned by profile, by gene
neighbourhood, or by system membership, and it cuts across the tree. `contextshift` takes
the labelling as an input and records where it came from, so a stale or disputed
classification scheme is visible rather than silent.

## The result is a 2×2, not a list of p-values

Divergence alone cannot say whether a site matters. Conservation alone cannot see a
between-group shift, because it averages over the groups. The cross of the two is the
result:

|                                      | not divergent | divergent |
|--------------------------------------|---------------|-----------|
| conserved everywhere                 | `core` — catalytic; the negative control | `suspect` — usually misalignment |
| conserved within groups, differs between | —         | **`determinant`** — the answer |
| constrained in one group only        | —             | `relaxed` — Type-I |
| variable everywhere                  | `variable`    | `variable` |

A column whose inputs are missing is emitted with `status` set. Nothing is dropped
quietly, and low-occupancy columns are flagged rather than filtered.

## Install

```bash
uv venv --python 3.12
uv pip install -e ".[dev]"
contextshift doctor          # what is on PATH, and whether diverge4 imports
```

External tools are optional and checked at call time: MMseqs2, MAFFT, FoldMason,
IQ-TREE 2, Rate4Site, MEME.

Functional divergence needs [DIVERGE v4](https://github.com/zjupgx/diverge4)
(Gu et al., *MBE* 2025; MIT). Installing it is not straightforward, and the
following was established by building and running 4.1.0, not by reading its docs:

- the PyPI name is **`diverge`**, not `diverge4`
- the sdist is broken — `setup.py` reads a `requirements.txt` it does not ship
- the default `src/` tree is Windows-targeted; build from a clone using `src_linux/`
- on macOS that tree still needs three patches: `isfinite` under clang,
  `#define version` colliding with a pybind11 member, and glibc-only `<values.h>`
- `Gu99`, `Rvs` and `TypeOneAnalysis` raise `NameError` on construction — they
  call an undefined `get_colnames`. **All Type-I entry points are therefore
  unusable in 4.1.0.** `Type2` works.
- `Gu99Batch` is defined but not exported from the package

## Shape of a run

```
members  ──derep──▶ weights ──align──▶ MSA ──tree──▶ per-group Newick
                                        │                    │
                                    conserve             diverge
                                        └────── classify ────┘
                                                  │
                                          sites × class
```

Every stage boundary is a declared table (`contextshift schemas`), so you can enter
anywhere with your own data and be told immediately if it does not fit.

## Check power before anything expensive

```bash
contextshift make-partition labels.tsv context.json \
  --name context --source mytool-1.2 --scheme "some-classification-2025"
contextshift power context.json
```

Dereplication keeps cluster size as a weight, and group size is reported as Kish
effective N. Raw counts overstate independence whenever a clade has been sequenced
many times; a group that looks large and is not gets flagged as underpowered before
you spend anything on it.

## Design rules

- The library is target-agnostic. Domain terms in library source fail a test.
- Dropped-in reference data never substitutes for computed results — it produces a
  comparison, so "we agree with the prior study" cannot be confused with "we overwrote
  our answer with theirs".
- Multiple-testing correction is global by default, across family × pair × column.
- Skipped comparisons are returned, never silently omitted.

## What DIVERGE actually returns

Verified by running 4.1.0 on its own CASP test data:

- `.summary` and `.results` are **properties**, not methods
- tree files are **positional** arguments; the `trees=` keyword takes Bio.Phylo
  objects, not paths
- `.results` carries the alignment position in its **index** (named `Position`,
  0-based), not a column, and the index is **sparse** — only positions DIVERGE
  kept appear (781 of 2088 columns on CASP, indexed 132..1460)
- the value is a **posterior probability Qk in [0, 1]**, not a p-value. It must
  never be routed through a p-value correction. `theta` and `alpha` are
  per-comparison and live in `.summary`
- the results column is named by joining the cluster names with `/` (`"A/B"`)
- tree depth must be **strictly greater than 3**, measured as
  `max(len(tree.trace(root, leaf)))` — not nesting depth. The library's own
  error message says "less than 3", which is wrong. Our check is tested for
  agreement with `diverge.binding.check_tree` on both sides of the boundary.

## Status

Early. Partition algebra, schemas, the 2×2 join, BH correction, DIVERGE tree
conformance and the drop zone are implemented and tested (60 tests). The
`type`, `families`, `map`, `neighbours` and `report` stages are not yet built,
so `members`, `neighbours` and `mapping` currently have no producing code.

Parsers for Rate4Site and MEME output are written against documented formats
and are **not yet checked against real tool output**.

## License

MIT
