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
| MMseqs2 `_cluster.tsv` is representative-then-member, tab separated | documented | MMseqs2 wiki |
| MMseqs2 `--min-seq-id`, `-c`, `--cov-mode 0` | documented | MMseqs2 wiki |
| FoldMason `easy-msa` writes `<prefix>_aa.fa` | documented | FoldMason README, verbatim |
| Rate4Site is ConSurf's engine | documented | ConSurf 2010, *NAR* 38:W529; Pupko et al. 2002 |
| FlaGs `_operon.tsv` column order | measured | writer in `FlaGs.py`, plus its shipped example output |
| FlaGs `outdesc` is `cluster(count)\taccession\tdescription` | measured | shipped example output |
| MAFFT `--localpair --maxiterate 1000` | unchecked | not installed |
| IQ-TREE `-B` (UFBoot), `--prefix`, `-m MFP` | unchecked | `.treefile` is documented; the three flags are not verified from a primary source |
| MEME `meme.txt` MOTIF line layout | unchecked | regex written against the documented format, not against real output |
| Rate4Site `.res` layout | unchecked | parser written against the documented format, not against real output |

Anything marked unchecked should be confirmed against real tool output before
its numbers are trusted.
