"""AlphaFold model retrieval and its two gates: sequence identity and confidence."""

import urllib.error
import urllib.request

import pytest

from contextshift.stages import structures


def online() -> bool:
    try:
        urllib.request.urlopen("https://alphafold.ebi.ac.uk/api/prediction/Q97ZJ4", timeout=15)
        return True
    except Exception:
        return False


needs_net = pytest.mark.skipif(not online(), reason="AlphaFold DB unreachable")


def test_confidence_gate():
    low = structures.Model("X", 6, 42.0, "AAA", "")
    high = structures.Model("X", 6, 92.0, "AAA", "")
    assert not low.confident and high.confident
    assert structures.check(low, "AAA") == ["mean pLDDT 42.0 < 70.0"]
    assert structures.check(high, "AAA") == []


def test_sequence_mismatch_is_reported_by_kind():
    m = structures.Model("X", 6, 95.0, "MKVL", "")
    assert "length" in structures.check(m, "MK")[0]
    assert "differs at 2 positions" in structures.check(m, "MKQQ")[0]


def test_no_expected_sequence_skips_that_check():
    m = structures.Model("X", 6, 95.0, "MKVL", "")
    assert structures.check(m, "") == []


@needs_net
def test_lookup_returns_a_versioned_model():
    m = structures.lookup("Q97ZJ4")
    assert m is not None
    # v3 and v4 now 404; the version must come from the API, not be constructed
    assert m.version >= 6
    assert m.cif_url.endswith(f"model_v{m.version}.cif")
    assert m.sequence and 0.0 < m.plddt <= 100.0


@needs_net
@pytest.mark.parametrize("acc", ["NOTAREALACC123", "ZZZZZZ"])
def test_missing_or_malformed_accessions_return_none(acc):
    # AFDB answers 404 for absent and 400 for malformed; neither should raise
    assert structures.lookup(acc) is None


@needs_net
def test_fetch_writes_a_cif_and_caches(tmp_path):
    m = structures.fetch("Q97ZJ4", tmp_path)
    assert m is not None and m.path.exists()
    assert m.path.read_text().lstrip().startswith("data_")
    first = m.path.stat().st_mtime_ns
    again = structures.fetch("Q97ZJ4", tmp_path)
    assert again.path.stat().st_mtime_ns == first


@needs_net
def test_batch_reports_every_miss(tmp_path):
    df, notes = structures.fetch_many(["Q97ZJ4", "NOTAREALACC123"], tmp_path)
    assert list(df["accession"]) == ["Q97ZJ4"]
    assert any("no AlphaFold model" in n for n in notes)
