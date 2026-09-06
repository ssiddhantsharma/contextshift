# Method provenance

What each stage does, the method it implements, and the source for that method.
Choices marked *substituted* differ from a tool named in the original plan;
the reason is given.

| Stage | Method | Source |
|---|---|---|
| `loci` | Detect and subtype systems, then explode loci into per-gene hits | CRISPRCasTyper: Russell et al. 2020, *CRISPR J* 3:462 |
| `loci.gene_presence` | Measure gene presence per subtype rather than assume it | — (guards against inheriting a review figure) |
| `families` | Split fused ORFs; index paralogous copies within a locus | — (mechanical; the failure modes are documented in tests) |
| `derep` | Cluster to remove pseudo-replication, keeping cluster size as a weight | MMseqs2: Steinegger & Söding 2017, *Nat Biotechnol* 35:1026 |
| `partition.effective_counts` | Effective sample size as summed weights (clusters, not sequences) | — (Kish effective N rejected: with equal weights it returns the raw count) |
| `align.mafft` | Progressive/iterative MSA | MAFFT: Katoh & Standley 2013, *MBE* 30:772 |
| `align.foldmason` | Structure-guided MSA over 3Di | FoldMason: Gilchrist et al. 2024, preprint; Foldseek 3Di: van Kempen et al. 2024, *Nat Biotechnol* 42:243 |
| `tree` | ML phylogeny with model selection and UFBoot supports | IQ-TREE 2: Minh et al. 2020, *MBE* 37:1530; ModelFinder: Kalyaanamoorthy et al. 2017, *Nat Methods* 14:587; UFBoot2: Hoang et al. 2018, *MBE* 35:518 |
| `conserve` (default) | Per-column Jensen-Shannon divergence from a background distribution, with Henikoff sequence weighting | *substituted for Rate4Site.* Capra & Singh 2007, *Bioinformatics* 23:1875; weighting: Henikoff & Henikoff 1994, *JMB* 243:574 |
| `conserve.rate4site` | Model-based per-site evolutionary rate | Rate4Site: Pupko et al. 2002, *Bioinformatics* 18:S71; as used by ConSurf: Ashkenazy et al. 2010, *NAR* 38:W529 |
| `diverge` | Type-I and Type-II functional divergence between subfamilies | Gu 1999, *MBE* 16:1664; Gu 2006, *MBE* 23:1937; DIVERGE v4: Cheng et al. 2025, *MBE* 42:msaf277 |
| `stats` | Benjamini-Hochberg FDR, global across family x pair x column | Benjamini & Hochberg 1995, *JRSS-B* 57:289 |
| `join` | Cross conservation with divergence to separate shared function from specificity | — (the 2x2 is this library's contribution) |
| `mapping` | Alignment column to author residue number, with distances to ligands and other chains | biotite: Kunzmann & Hamacher 2018, *BMC Bioinformatics* 19:346 |
| `motifs` | Ungapped motif discovery per group | MEME: Bailey & Elkan 1994, *ISMB* 2:28 |
| `neighbours` | Gene neighbourhood conservation | FlaGs: Saha et al. 2021, *Bioinformatics* 37:1312 |
| `reference` | Score the pipeline against a prior study's published claims | — |

## Substitutions

**Rate4Site to Jensen-Shannon.** Rate4Site is not packaged for Homebrew and was
not run, so its output format would have stayed an unverified assumption. Capra
& Singh (2007) benchmarked conservation measures for identifying functional
sites and found JS divergence the strongest of the simple scores; it needs no
external tool, so it is unit-tested here. Rate4Site remains available for a
model-based rate. Both report on the same convention: lower is more constrained.

**MEME is optional.** On an alignable family, per-column conservation and the
divergence tests already carry the signal that motif discovery would restate.
MEME earns its place on sets that are hard to align, and its adapter is kept
for that, unrun.

## What this pipeline does not do

- It does not infer causality: a divergent site is a hypothesis for assay.
- It does not correct a bad genome set. If the input is not dereplicated at the
  genome level, no later stage can repair the resulting sampling bias.
- Type-I results depend on a shim for a defect in DIVERGE 4.1.0, and every run
  that used it says so.
