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
IQ-TREE 2, Rate4Site, MEME. Functional divergence needs
[diverge4](https://github.com/zjupgx/diverge4) (Gu et al., *MBE* 2025).

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

## Status

Early. The partition algebra, schemas, 2×2 join, FDR, tree conformance and drop-zone are
implemented and tested. `stages.diverge.normalise` is written against DIVERGE's documented
column names and **has not yet been validated against installed diverge4 output** — verify
it before trusting any number it produces.

## License

MIT
