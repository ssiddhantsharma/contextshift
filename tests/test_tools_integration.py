"""Adapters run against the real tools, where those tools are installed."""

from pathlib import Path

import pytest

from contextshift.stages import align, derep
from contextshift.stages._external import MAFFT, MMSEQS

CASP = Path(__file__).parent / "data" / "diverge_casp" / "CASP.aln"


@pytest.fixture(scope="module")
def unaligned(tmp_path_factory):
    from Bio import AlignIO

    d = tmp_path_factory.mktemp("seqs")
    p = d / "seqs.fa"
    aln = AlignIO.read(str(CASP), "fasta")
    with open(p, "w") as f:
        for rec in aln[:8]:
            f.write(f">{rec.id}\n{str(rec.seq).replace('-', '')[:120]}\n")
    return p


@pytest.mark.skipif(MAFFT.path is None, reason="mafft not installed")
def test_mafft_invocation_produces_an_alignment(unaligned, tmp_path):
    out = align.mafft(unaligned, tmp_path / "aln.fa", threads=2)
    occ = align.occupancy(out)
    assert out.exists()
    assert len(occ) > 0
    assert 0.0 < min(occ) <= max(occ) <= 1.0


@pytest.mark.skipif(MMSEQS.path is None, reason="mmseqs not installed")
def test_mmseqs_cluster_tsv_is_representative_then_member(unaligned, tmp_path):
    d = derep.cluster(unaligned, identity=0.3, threads=2, workdir=tmp_path / "mm")
    assert set(d.columns) >= {"member_id", "representative_id", "cluster_size", "weight"}
    # every representative is also a member of its own cluster
    assert set(d["representative_id"]).issubset(set(d["member_id"]))


@pytest.mark.skipif(MMSEQS.path is None, reason="mmseqs not installed")
def test_weights_sum_to_the_number_of_clusters(unaligned, tmp_path):
    d = derep.cluster(unaligned, identity=0.3, threads=2, workdir=tmp_path / "mm2")
    assert float(d["weight"].sum()) == pytest.approx(d["cluster_id"].nunique())


@pytest.mark.skipif(MMSEQS.path is None, reason="mmseqs not installed")
def test_effective_n_counts_clusters_not_sequences(unaligned, tmp_path):
    from contextshift.partition import Partition, Provenance

    d = derep.cluster(unaligned, identity=0.3, threads=2, workdir=tmp_path / "mm3")
    p = Partition(
        name="x",
        labels={m: "all" for m in d["member_id"]},
        provenance=Provenance(source="test"),
        weights=derep.member_weights(d),
    )
    assert p.counts()["all"] == len(d)
    assert p.effective_counts()["all"] == pytest.approx(d["cluster_id"].nunique())
