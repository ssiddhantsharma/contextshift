"""A labelling of family members into groups to be contrasted.

Labels come from outside the phylogeny, never from cutting the tree the
analysis is about to test.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from itertools import combinations
from math import comb
from pathlib import Path


@dataclass(frozen=True)
class Provenance:
    """Where labels came from."""

    source: str
    version: str = "unknown"
    scheme: str = "unknown"
    notes: str = ""

    def as_dict(self) -> dict[str, str]:
        return {
            "source": self.source,
            "version": self.version,
            "scheme": self.scheme,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class GroupPower:
    label: str
    n: int
    effective_n: float
    caveat: str = ""

    @property
    def underpowered(self) -> bool:
        return self.effective_n < MIN_EFFECTIVE_N

    @property
    def qualified(self) -> bool:
        """True if anything about this group needs stating alongside its result."""
        return self.underpowered or bool(self.caveat)


MIN_EFFECTIVE_N = 30.0


@dataclass
class Partition:
    name: str
    labels: dict[str, str]
    provenance: Provenance
    weights: dict[str, float] = field(default_factory=dict)
    #: Per-group qualifications carried through to reports.
    caveats: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.labels:
            raise ValueError(f"partition {self.name!r} has no labels")
        unknown = set(self.weights) - set(self.labels)
        if unknown:
            raise ValueError(f"weights reference {len(unknown)} members not in labels")
        stray = set(self.caveats) - set(self.labels.values())
        if stray:
            raise ValueError(f"caveats reference unknown groups: {sorted(stray)}")

    @property
    def members(self) -> list[str]:
        return list(self.labels)

    @property
    def group_names(self) -> list[str]:
        return sorted(set(self.labels.values()))

    def weight(self, member: str) -> float:
        return self.weights.get(member, 1.0)

    def groups(self) -> dict[str, list[str]]:
        out: dict[str, list[str]] = {g: [] for g in self.group_names}
        for member, label in self.labels.items():
            out[label].append(member)
        return {g: sorted(m) for g, m in out.items()}

    def counts(self) -> dict[str, int]:
        return dict(Counter(self.labels.values()))

    def effective_counts(self) -> dict[str, float]:
        """Independent evidence per group: summed weights.

        Dereplication gives each cluster one unit shared among its members, so
        this counts clusters, not sequences. Kish effective N is not used: with
        equal weights it returns the raw count.
        """
        return {
            label: sum(self.weight(m) for m in members)
            for label, members in self.groups().items()
        }

    def power(self) -> list[GroupPower]:
        eff = self.effective_counts()
        counts = self.counts()
        return [
            GroupPower(g, counts[g], eff[g], self.caveats.get(g, ""))
            for g in self.group_names
        ]

    def qualified_groups(self) -> dict[str, str]:
        """Groups carrying a caveat."""
        return {g: c for g, c in self.caveats.items() if c}

    def underpowered_groups(self) -> list[str]:
        return [p.label for p in self.power() if p.underpowered]

    def pairs(self) -> list[tuple[str, str]]:
        return list(combinations(self.group_names, 2))

    def restrict(self, members: list[str]) -> Partition:
        keep = set(members) & set(self.labels)
        return Partition(
            name=self.name,
            labels={m: self.labels[m] for m in sorted(keep)},
            provenance=self.provenance,
            weights={m: w for m, w in self.weights.items() if m in keep},
            caveats=dict(self.caveats),
        )

    def drop_groups(self, labels: list[str]) -> Partition:
        drop = set(labels)
        return self.restrict([m for m, g in self.labels.items() if g not in drop])

    def cross(self, other: Partition, sep: str = "|") -> Partition:
        """Cell-wise product over shared members."""
        shared = sorted(set(self.labels) & set(other.labels))
        if not shared:
            raise ValueError(f"{self.name} and {other.name} share no members")
        return Partition(
            name=f"{self.name}{sep}{other.name}",
            labels={m: f"{self.labels[m]}{sep}{other.labels[m]}" for m in shared},
            provenance=Provenance(
                source=f"cross({self.provenance.source},{other.provenance.source})",
                scheme=f"{self.provenance.scheme}{sep}{other.provenance.scheme}",
            ),
            weights={m: self.weight(m) for m in shared},
        )

    def to_json(self, path: Path) -> None:
        path.write_text(
            json.dumps(
                {
                    "name": self.name,
                    "provenance": self.provenance.as_dict(),
                    "labels": self.labels,
                    "weights": self.weights,
                    "caveats": self.caveats,
                },
                indent=2,
                sort_keys=True,
            )
        )

    @classmethod
    def from_json(cls, path: Path) -> Partition:
        d = json.loads(Path(path).read_text())
        return cls(
            name=d["name"],
            labels=d["labels"],
            provenance=Provenance(**d["provenance"]),
            weights=d.get("weights", {}),
            caveats=d.get("caveats", {}),
        )

    @classmethod
    def from_tsv(cls, path: Path, name: str, provenance: Provenance) -> Partition:
        labels: dict[str, str] = {}
        weights: dict[str, float] = {}
        for line in Path(path).read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            labels[fields[0]] = fields[1]
            if len(fields) > 2:
                weights[fields[0]] = float(fields[2])
        return cls(name=name, labels=labels, provenance=provenance, weights=weights)


def adjusted_rand_index(a: Partition, b: Partition) -> float:
    """ARI over shared members; 1.0 is identical grouping."""
    shared = sorted(set(a.labels) & set(b.labels))
    n = len(shared)
    if n < 2:
        return float("nan")

    table: Counter[tuple[str, str]] = Counter(
        (a.labels[m], b.labels[m]) for m in shared
    )
    row: Counter[str] = Counter()
    col: Counter[str] = Counter()
    for (ra, cb), v in table.items():
        row[ra] += v
        col[cb] += v

    sum_ij = sum(comb(v, 2) for v in table.values())
    sum_i = sum(comb(v, 2) for v in row.values())
    sum_j = sum(comb(v, 2) for v in col.values())
    total = comb(n, 2)

    expected = sum_i * sum_j / total
    maximum = 0.5 * (sum_i + sum_j)
    if maximum == expected:
        return 1.0
    return (sum_ij - expected) / (maximum - expected)
