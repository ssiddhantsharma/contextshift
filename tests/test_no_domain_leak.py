"""The library stays target-agnostic."""

import re
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src" / "contextshift"

# Biological findings, not tool names. Naming a tool the library drives is
# fine: mafft, mmseqs and cctyper are dependencies, cas4 is a research subject.
FORBIDDEN = [
    "cas4", "cas9", "cas12", "csa1", "crispr", "solo-cas", "recb",
    "protospacer", "spacer acquisition", "hudaiberdiev", "makarova",
]


def library_files():
    return sorted(p for p in SRC.rglob("*.py"))


def test_no_domain_terms_in_library():
    offences = []
    for path in library_files():
        rel = path.relative_to(SRC).as_posix()
        text = path.read_text().lower()
        for term in FORBIDDEN:
            for m in re.finditer(re.escape(term), text):
                line = text[: m.start()].count("\n") + 1
                offences.append(f"{rel}:{line} contains {term!r}")
    assert not offences, "domain terms leaked into the library:\n" + "\n".join(offences)


def test_library_has_no_hardcoded_paths():
    offences = []
    for path in library_files():
        for i, line in enumerate(path.read_text().splitlines(), start=1):
            if "/Users/" in line or "/home/" in line:
                offences.append(f"{path.relative_to(SRC)}:{i}")
    assert not offences, f"absolute paths in library: {offences}"
