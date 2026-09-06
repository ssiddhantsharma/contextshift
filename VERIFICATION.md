# Verification status of external contracts

Every assumption this library makes about another tool's interface or output,
and how it was checked. "Measured" means the tool was run or its source read;
"documented" means a primary source states it; "unchecked" means neither.

| Contract | Status | Evidence |
|---|---|---|
| DIVERGE tree files are positional; `trees=` takes Bio.Phylo objects | measured | `binding.py` `BaseAnalysis.__init__` |
| DIVERGE `.summary` / `.results` are properties | measured | ran 4.1.0 |
| DIVERGE position is the results index, sparse | measured | CASP: 781 of 2088 columns, 132..1460 |
| DIVERGE value is posterior Qk in [0, 1] | measured | CASP: 0.0–0.849 |
| DIVERGE result column is `"A/B"` | measured | `_r_names()` returns `['A/B']` |
| DIVERGE tree depth must be `> 3`, counted as root-to-leaf edges | measured | agrees with `binding.check_tree` on both sides of the boundary |
| `Gu99`/`Rvs`/`TypeOneAnalysis` raise `NameError` in 4.1.0 | measured | `get_colnames` undefined; upstream PR zjupgx/diverge4#8 |
| Qk ≥ 0.9 as the calling threshold | documented | DIVERGE User Guide uses `results.iloc[:, 0] > 0.9` |
| MMseqs2 `_cluster.tsv` is representative-then-member, tab separated | measured | ran 18-8cc5c; every representative is also a member |
| Dereplication weights sum to the cluster count | measured | asserted in `test_tools_integration.py` |
| MMseqs2 `--min-seq-id`, `-c`, `--cov-mode 0` | documented | MMseqs2 wiki |
| FoldMason `easy-msa` writes `<prefix>_aa.fa` | documented | FoldMason README, verbatim |
| Rate4Site is ConSurf's engine | documented | ConSurf 2010, *NAR* 38:W529; Pupko et al. 2002 |
| FlaGs `_operon.tsv` column order | measured | writer in `FlaGs.py`, plus its shipped example output |
| FlaGs `outdesc` is `cluster(count)\taccession\tdescription` | measured | shipped example output |
| CCTyper `cas_operons.tab` list columns are stringified Python lists | measured | `castyping.py` writes `list(tmp['Hmm'])` into a cell, then `to_csv` |
| CCTyper column names | documented | its README, cross-checked against the dict in `castyping.py` |
| MAFFT `--localpair --maxiterate 1000 --anysymbol` | measured | ran 7.526 through the adapter |
| IQ-TREE `-B` (UFBoot, `>=1000`), `--prefix`, `-m MFP` | measured | `-h` on 2.4.0 (macOS) and 3.1.3 (container); full run wrote `.treefile` with UFBoot supports |
| IQ-TREE binary name varies (`iqtree2`/`iqtree3`/`iqtree`) | measured | bioconda ships `iqtree3`; `Tool.aliases` resolves it |
| MEME `meme.txt` MOTIF line layout | unchecked, optional | not packaged for Homebrew; not run. Motif discovery is optional, and per-group comparison is also obtainable from the per-scope conservation |
| Jensen-Shannon conservation | measured | computed in-library; unit-tested for the conserved/variable ordering, gap handling and redundancy weighting |
| Rate4Site `.res` layout | unchecked, optional | parser retained as an alternative; not run. The default conservation path no longer needs it |

Anything marked unchecked should be confirmed against real tool output before
its numbers are trusted.

## Container

`Dockerfile` builds DIVERGE from source and installs MMseqs2, MAFFT, IQ-TREE
and MEME. Verified by building and running it: DIVERGE 4.1.0 compiles and
imports, `needs_shim` is False, and `doctor` finds four of the optional tools.

Three bugs were found only by building it, not by reading it: heredocs are not
honoured inside `RUN`; `WORKDIR` creates directories as root so the unprivileged
user cannot write to them; and `MAMBA_DOCKERFILE_ACTIVATE` applies at build time
only, so the runtime entrypoint must activate the environment.

CCTyper is not in the image on arm64. The cause is not a missing ARM build:
prodigal 2.6.3, minced, cairosvg, blast and hmmer are all available for
linux-aarch64 or noarch. The bioconda cctyper recipe pins
`prodigal >=2.0,<=2.6.2`, and 2.6.2 is the last version *without* an aarch64
build, so the ceiling excludes the only ARM-capable prodigal. Verified by
resolving the same dependency set with the ceiling removed, which succeeds.
Fix proposed as bioconda-recipes#68871.
