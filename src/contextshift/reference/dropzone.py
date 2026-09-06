"""Reference artifacts from a prior study, used to score this pipeline.

Dropped files never substitute for computed results; each slot produces a
comparison.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..schema import TableSchema

Comparator = Callable[[Any, Any], "Comparison"]


@dataclass(frozen=True)
class Comparison:
    slot: str
    metric: str
    value: float
    expected: float | None = None
    passed: bool | None = None
    detail: str = ""


@dataclass(frozen=True)
class Slot:
    name: str
    filename: str
    description: str
    schema: TableSchema | None = None
    loader: Callable[[Path], Any] | None = None
    comparator: Comparator | None = None

    def load(self, root: Path) -> Any:
        path = root / self.filename
        if not path.exists():
            return None
        if self.loader is not None:
            return self.loader(path)
        if self.schema is not None:
            return self.schema.read(path)
        return path.read_text()


@dataclass
class DropZone:
    root: Path
    slots: list[Slot] = field(default_factory=list)

    def register(self, slot: Slot) -> None:
        if any(s.name == slot.name for s in self.slots):
            raise ValueError(f"slot {slot.name!r} already registered")
        self.slots.append(slot)

    @property
    def expected_files(self) -> set[str]:
        return {s.filename for s in self.slots}

    def present(self) -> list[Slot]:
        return [s for s in self.slots if (self.root / s.filename).exists()]

    def absent(self) -> list[Slot]:
        return [s for s in self.slots if not (self.root / s.filename).exists()]

    def unmapped(self) -> list[Path]:
        """Files in the drop zone that no slot claims."""
        if not self.root.exists():
            return []
        return sorted(
            p
            for p in self.root.rglob("*")
            if p.is_file() and p.relative_to(self.root).as_posix() not in self.expected_files
        )

    def compare(self, computed: dict[str, Any]) -> list[Comparison]:
        out: list[Comparison] = []
        for slot in self.present():
            if slot.comparator is None or slot.name not in computed:
                continue
            reference = slot.load(self.root)
            out.append(slot.comparator(computed[slot.name], reference))
        return out

    def report(self, computed: dict[str, Any] | None = None) -> str:
        lines = [f"drop zone: {self.root}"]
        present, absent = self.present(), self.absent()
        lines.append(f"  slots filled : {len(present)}/{len(self.slots)}")
        for slot in present:
            lines.append(f"    [x] {slot.name:24s} {slot.filename}")
        for slot in absent:
            lines.append(f"    [ ] {slot.name:24s} {slot.filename}  ({slot.description})")

        extra = self.unmapped()
        if extra:
            lines.append(f"  unmapped files: {len(extra)}")
            for p in extra:
                lines.append(f"    ??  {p.relative_to(self.root)}")

        if computed:
            comparisons = self.compare(computed)
            lines.append(f"  comparisons  : {len(comparisons)}")
            for c in comparisons:
                mark = "" if c.passed is None else ("PASS" if c.passed else "FAIL")
                expected = "" if c.expected is None else f" (expected {c.expected:g})"
                lines.append(
                    f"    {mark:4s} {c.slot}.{c.metric} = {c.value:.4g}{expected} {c.detail}".rstrip()
                )
        return "\n".join(lines)
