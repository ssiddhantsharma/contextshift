# Test fixtures

## `diverge_casp/`
`CASP.aln`, `cl1.tree`, `cl2.tree` — from `Tutorial/test_data` of
[zjupgx/diverge4](https://github.com/zjupgx/diverge4) (MIT, © 2025 ZjuPgx).
Used so the DIVERGE adapter is tested against the inputs DIVERGE ships.

## `structures/4IC1_A.cif`
Chain A of [PDB 4IC1](https://www.rcsb.org/structure/4IC1), waters removed.
X-ray, 2.35 Å. Carries an SF4 cluster and a catalytic Mn, so residue anchors
in the mapping tests are measured from coordinates rather than assumed.

## `catalytic/`
Five alignments with Catalytic Site Atlas labels, from the data supporting
Capra & Singh 2007, *Bioinformatics* 23:1875-1882
(compbio.cs.princeton.edu/conservation/). Alignments converted from CLUSTAL to
FASTA; labels unchanged. Used to check the conservation stage against published
ground truth rather than against its own definition.
