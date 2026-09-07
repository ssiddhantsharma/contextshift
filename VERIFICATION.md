# Verification

What this library assumes about the outside world, and how each assumption was
checked. Kept because several assumptions taken from documentation turned out
to be wrong, and each would have corrupted results silently.

**measured** — the tool was run, or its source read.
**documented** — a primary source states it.
**unchecked** — neither.

## Validated against published benchmarks

Data supporting Capra & Singh 2007, *Bioinformatics* 23:1875-1882
(compbio.cs.princeton.edu/conservation/), whose Jensen-Shannon method the
conservation stage implements. Reproduce with:

```bash
python benchmarks/capra_singh.py --workdir /tmp/cs
```

| Check | Result |
|---|---|
| Do catalytic sites score as more constrained? | **57 of 57** proteins; mean per-protein AUC **0.924** |
| Do residue-to-ligand distances match published values? | **7 of 8** structures exact (max difference **0.00 A**) |

The eighth, 1ADB, differs by up to 1.73 A and has **1,410 deposited hydrogens**;
the other structures have none. One residue is published at 1.59 A, shorter than
a covalent bond, so the published distances include hydrogens. This stage
measures heavy-atom distances deliberately, and the difference is about one
bond length.

Small fixtures from both ship in `tests/data/`, so the tests run offline.

Two things the ligand comparison established about calling `distance_to`:

- The nearest ligand atom for a residue in one chain is often a **copy of the
  ligand bound to another chain**. Restricting the target to a single chain is
  wrong: on 1ADB it moved the maximum deviation from 1.7 A to 33 A.
- Per-protein AUC is noise below about three scored catalytic sites, so the
  catalytic fixture requires that and the aggregate test pools positions.

Not usable as shipped: the benchmark's **interface labels** are indexed by
alignment position, and those alignments are not included (its README points to
Caffrey et al. 2004). 128 label files, no alignments.

## Internal consistency

| Contract | Status | Evidence |
|---|---|---|
| A pairwise divergence estimate is unaffected by the other clusters | measured | theta and every per-site posterior byte-identical between a 2-cluster run and the same pair from a 3-cluster run, over the same 781 kept positions. The library loops over pairs, so this had to hold |
| Dereplication weights sum to the cluster count | measured | asserted in `test_tools_integration.py` |

## DIVERGE

Built from source and run; none of this is from its documentation.

| Contract | Status | Evidence |
|---|---|---|
| Tree files are positional; `trees=` takes Bio.Phylo objects, not paths | measured | `binding.py` `BaseAnalysis.__init__` |
| `.summary` and `.results` are properties, not methods | measured | ran 4.1.0 |
| Position is the results **index**, and is **sparse** | measured | CASP: 781 of 2088 columns, indexed 132..1460 |
| The value is a posterior **Qk in [0, 1]**, not a p-value | measured | CASP range 0.0-0.849. It must never be fed to a p-value correction |
| Result column is named `"A/B"` from the cluster names | measured | `_r_names()` returns `['A/B']` |
| Tree depth must be **> 3**, counted as root-to-leaf edges | measured | agrees with `binding.check_tree` either side of the boundary. Not nesting depth, and the library's own error message says "less than 3" |
| `Gu99`, `Rvs`, `TypeOneAnalysis` raise `NameError` in 4.1.0 | measured | undefined `get_colnames`; every Type-I entry point is unusable as shipped. Shim recorded in run notes. Upstream [PR #8](https://github.com/zjupgx/diverge4/pull/8) |
| Qk >= 0.9 as the calling threshold | documented | DIVERGE User Guide |
| Not installable from PyPI | measured | sdist reads a `requirements.txt` it does not ship; default `src/` tree is MSVC-only. PRs [#9](https://github.com/zjupgx/diverge4/pull/9), [#10](https://github.com/zjupgx/diverge4/pull/10) |

## Other tools

| Contract | Status | Evidence |
|---|---|---|
| MMseqs2 `_cluster.tsv` is representative-then-member | measured | ran 18-8cc5c; every representative is also a member |
| MMseqs2 `--min-seq-id`, `-c`, `--cov-mode 0` | documented | MMseqs2 wiki |
| MAFFT `--localpair --maxiterate 1000 --anysymbol` | measured | ran 7.526 through the adapter |
| IQ-TREE `-B` (UFBoot, `>=1000`), `--prefix`, `-m MFP` | measured | `-h` on 2.4.0 and 3.1.3; a run wrote `.treefile` with UFBoot supports |
| IQ-TREE binary name varies (`iqtree2`/`iqtree3`/`iqtree`) | measured | bioconda ships `iqtree3`; `Tool.aliases` resolves it |
| FoldMason `easy-msa` writes `<prefix>_aa.fa` | documented | README, verbatim |
| Rate4Site is ConSurf's engine | documented | ConSurf 2010, *NAR* 38:W529; Pupko et al. 2002 |
| CCTyper `cas_operons.tab` list columns are stringified Python lists | measured | `castyping.py` writes `list(tmp['Hmm'])` into a cell, then `to_csv` |
| CCTyper column names | documented | its README, cross-checked against the dict in `castyping.py` |
| DefenseFinder subtype string is `CAS_Class1-Subtype-I-E` | documented | raw README, verbatim |
| DefenseFinder `protein_in_syst` and `name_of_profiles_in_sys` are independently alphabetised and must **not** be zipped | documented | raw README warning, verbatim |
| PADLOC systems are `cas_type_I-E`, and it subdivides I-B1/I-B2, I-F1/2/3 | measured | padloc-db `sys/` listing, 33 CRISPR systems |
| FlaGs `_operon.tsv` column order | measured | the writer in `FlaGs.py`, plus its shipped example |
| pyhmmer returns `str`, not `bytes`, from 0.11 onward | measured | ran 0.12.3 |
| Pfam HMMs come gzipped from InterPro; pfam.xfam.org is retired | measured | PF01930 = `Cas_Cas4`, PF06023 = `Csa1` |
| AlphaFold DB is on model v6; v3 and v4 return 404 | measured | version resolved via the API, never constructed |
| AlphaFold DB answers 400 for malformed and 404 for absent accessions | measured | both treated as "no model" |
| SpacePHARER writes `#` match lines and `>` hit lines, 9 fields, PAM pipe-separated | documented | raw README, verbatim |
| MEME `meme.txt` MOTIF line layout | unchecked, optional | not packaged for Homebrew; not run |
| Rate4Site `.res` layout | unchecked, optional | parser kept as an alternative; the default conservation path does not need it |

## Container

`Dockerfile` builds DIVERGE from source and installs MMseqs2, MAFFT, IQ-TREE,
MEME, pyhmmer and DefenseFinder. Verified by building and running it: DIVERGE
4.1.0 compiles and imports, `needs_shim` is False, pyhmmer 0.12.3 imports and
`defense-finder` is on PATH.

Four bugs were found only by building it: heredocs are not honoured inside
`RUN`; `WORKDIR` creates directories as root so an unprivileged user cannot
write to them; `MAMBA_DOCKERFILE_ACTIVATE` applies at build time only, so the
runtime entrypoint must activate the environment; and bioconda ships `iqtree3`,
which the tool detection did not know about.

PADLOC is left out of the default image. It resolves on arm64, but its R stack
takes the install from 154 packages / 341MB to 330 / 695MB. One line adds it.

CCTyper is absent on arm64, and not because a dependency lacks an ARM build --
prodigal 2.6.3, minced, cairosvg, blast and hmmer are all available for
linux-aarch64 or noarch. The bioconda recipe pins `prodigal >=2.0,<=2.6.2`, and
2.6.2 is the last version *without* an aarch64 build, so the ceiling excludes
the only ARM-capable prodigal. Verified by re-resolving with the ceiling
removed, which succeeds. Fix proposed as
[bioconda-recipes #68871](https://github.com/bioconda/bioconda-recipes/pull/68871).
