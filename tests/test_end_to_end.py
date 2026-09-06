"""The stages chain: typed loci in, classified sites with structural context out."""

from pathlib import Path

import pandas as pd
import pytest

from contextshift import report
from contextshift.join import classify
from contextshift.partition import Partition, Provenance
from contextshift.stages import conserve, derep, diverge, families, loci, mapping
from contextshift.stages._external import MMSEQS

DATA = Path(__file__).parent / "data"
CASP = DATA / "diverge_casp"
CIF = DATA / "structures" / "4IC1_A.cif"

pytestmark = pytest.mark.skipif(not diverge.available(), reason="diverge not installed")


@pytest.fixture(scope="module")
def operons(tmp_path_factory):
    """Two loci per subtype, one gene each, mirroring cas_operons.tab."""
    p = tmp_path_factory.mktemp("loci") / "cas_operons.tab"
    rows = ["Contig\tOperon\tStart\tEnd\tPrediction\tGenes\tPositions"]
    for i, sub in enumerate(["I-A", "I-A", "I-C", "I-C"]):
        rows.append(
            f"c{i}\tc{i}@1\t1\t900\t{sub}\t{['GeneX', 'GeneY']}\t{[f'p{i}a', f'p{i}b']}"
        )
    p.write_text("\n".join(rows) + "\n")
    return p


def test_loci_to_partition_to_power(operons):
    hits = loci.to_hits(loci.read_operons(operons))
    members, freport = families.build(hits)
    partition, dropped = loci.subtype_partition(hits, "typer-1", "scheme-2020")

    assert freport.n_members == 8
    assert set(partition.group_names) == {"I-A", "I-C"}
    assert dropped == []
    # power is computable straight off the partition
    assert {r.label for r in partition.power()} == {"I-A", "I-C"}


@pytest.mark.skipif(MMSEQS.path is None, reason="mmseqs not installed")
def test_derep_weights_feed_effective_counts(tmp_path):
    from Bio import AlignIO

    aln = AlignIO.read(str(CASP / "CASP.aln"), "fasta")
    fa = tmp_path / "seqs.fa"
    with open(fa, "w") as f:
        for rec in aln[:10]:
            f.write(f">{rec.id}\n{str(rec.seq).replace('-', '')[:150]}\n")

    d = derep.cluster(fa, identity=0.3, threads=2, workdir=tmp_path / "mm")
    p = Partition(
        name="context",
        labels={m: "all" for m in d["member_id"]},
        provenance=Provenance(source="test"),
        weights=derep.member_weights(d),
    )
    assert p.effective_counts()["all"] <= p.counts()["all"]


def test_diverge_to_conservation_to_classified(tmp_path):
    """The analytical core: divergence x conservation on the same columns."""
    sites, comparisons, notes = diverge.run_pair(
        CASP / "CASP.aln", CASP / "cl1.tree", CASP / "cl2.tree",
        family="CASP", partition="context", group_a="A", group_b="B",
    )
    cons = conserve.jensen_shannon(CASP / "CASP.aln", "CASP", conserve.SCOPE_ALL)
    for group in ("A", "B"):
        cons = pd.concat(
            [cons, conserve.jensen_shannon(CASP / "CASP.aln", "CASP", group)],
            ignore_index=True,
        )

    out = classify(sites, cons)
    assert len(out) > 0
    # every classified column carries a class and a status, never a bare NaN
    assert out["class"].notna().all()
    assert out["status"].notna().all()
    # the columns DIVERGE reported are the ones classified
    assert set(out["column"]) <= set(sites["column"])


def test_mapping_attaches_structure_to_classified_columns():
    residues = mapping.read_chain(CIF, "A")
    seq = residues.sequence
    m, mismatches = mapping.map_columns(seq, residues, "Cas4like", "4IC1_A")
    d = mapping.distance_to(CIF, "A", residues, target_res_names=("SF4",))

    sites = pd.DataFrame([{
        "family": "Cas4like", "partition": "context", "group_a": "A", "group_b": "B",
        "column": int(m["column"].iloc[0]), "test": "type2", "posterior": 0.97,
    }])
    cons = pd.DataFrame([
        {"family": "Cas4like", "column": int(m["column"].iloc[0]), "scope": s, "rate": r}
        for s, r in [("all", 1.2), ("A", 0.1), ("B", 0.1)]
    ])
    annotated = mapping.annotate(classify(sites, cons), m, d)
    assert not mismatches
    assert annotated["resnum"].notna().all()
    assert annotated["min_distance"].notna().all()


def test_report_writes_from_a_classified_frame(tmp_path):
    sites, _, _ = diverge.run_pair(
        CASP / "CASP.aln", CASP / "cl1.tree", CASP / "cl2.tree",
        family="CASP", partition="context", group_a="A", group_b="B",
    )
    cons = conserve.jensen_shannon(CASP / "CASP.aln", "CASP", conserve.SCOPE_ALL)
    for g in ("A", "B"):
        cons = pd.concat([cons, conserve.jensen_shannon(CASP / "CASP.aln", "CASP", g)],
                         ignore_index=True)
    classified = classify(sites, cons)

    run = report.RunSummary()
    run.stage("sites", len(classified))
    written = report.write(tmp_path, run, classified)
    assert {p.name for p in written} == {"funnel.png", "classes.png", "summary.txt"}
    assert (tmp_path / "summary.txt").read_text().strip()
