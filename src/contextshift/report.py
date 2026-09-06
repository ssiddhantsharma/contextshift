"""Run summary: a funnel, per-group power, the class breakdown, and sites in context.

Every panel reports counts that were dropped or unresolved alongside those that
were kept, so a figure cannot hide an exclusion.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .join import CORE, DETERMINANT, RELAXED, SUSPECT, UNKNOWN, VARIABLE, summary
from .partition import Partition

CLASS_ORDER = [DETERMINANT, RELAXED, CORE, VARIABLE, SUSPECT, UNKNOWN]
PALETTE = {
    DETERMINANT: "#B03A2E",
    RELAXED: "#CA8A2E",
    CORE: "#2E6DB0",
    VARIABLE: "#B8BCC0",
    SUSPECT: "#6B3FA0",
    UNKNOWN: "#4D4D4D",
}


@dataclass
class RunSummary:
    stages: dict[str, int] = field(default_factory=dict)
    dropped: dict[str, int] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def stage(self, name: str, kept: int, dropped: int = 0) -> None:
        self.stages[name] = kept
        if dropped:
            self.dropped[name] = dropped

    def as_text(self) -> str:
        width = max((len(k) for k in self.stages), default=0)
        lines = ["stage counts"]
        for name, kept in self.stages.items():
            lost = self.dropped.get(name, 0)
            suffix = f"   (-{lost} dropped)" if lost else ""
            lines.append(f"  {name:<{width}}  {kept:>8d}{suffix}")
        for note in self.notes:
            lines.append(f"  note: {note}")
        return "\n".join(lines)


def _figure(nrows=1, ncols=1, size=(7, 4)):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt.subplots(nrows, ncols, figsize=size, constrained_layout=True)


def funnel(run: RunSummary, path: Path) -> Path:
    fig, ax = _figure(size=(7, 0.5 * max(len(run.stages), 3) + 1.5))
    names = list(run.stages)
    kept = [run.stages[n] for n in names]
    lost = [run.dropped.get(n, 0) for n in names]
    y = range(len(names))
    ax.barh(list(y), kept, color="#2E6DB0", label="kept")
    ax.barh(list(y), lost, left=kept, color="#B8BCC0", label="dropped")
    ax.set_yticks(list(y), names)
    ax.invert_yaxis()
    ax.set_xlabel("records")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("stage funnel")
    fig.savefig(path, dpi=150)
    return path


def power(partition: Partition, path: Path) -> Path:
    rows = partition.power()
    fig, ax = _figure(size=(7, 0.4 * max(len(rows), 3) + 1.5))
    y = range(len(rows))
    ax.barh(list(y), [r.n for r in rows], color="#B8BCC0", label="sequences")
    ax.barh(list(y), [r.effective_n for r in rows], color="#2E6DB0", label="effective (clusters)")
    labels = [f"{r.label}{' *' if r.caveat else ''}" for r in rows]
    ax.set_yticks(list(y), labels)
    ax.invert_yaxis()
    ax.set_xlabel("count")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title(f"power by {partition.name}   (* carries a caveat)")
    fig.savefig(path, dpi=150)
    return path


def classes(classified: pd.DataFrame, path: Path) -> Path:
    counts = summary(classified)
    pivot = (
        counts.groupby(["group_pair", "class"])["n"].sum().unstack(fill_value=0)
        if "group_pair" in counts.columns
        else counts.groupby("class")["n"].sum().to_frame().T
    )
    present = [c for c in CLASS_ORDER if c in pivot.columns]
    fig, ax = _figure(size=(7, 4))
    bottom = [0] * len(pivot)
    for cls in present:
        vals = list(pivot[cls])
        ax.bar(range(len(pivot)), vals, bottom=bottom, label=cls, color=PALETTE[cls])
        bottom = [b + v for b, v in zip(bottom, vals, strict=True)]
    ax.set_xticks(range(len(pivot)), [str(i) for i in pivot.index], rotation=30, ha="right")
    ax.set_ylabel("alignment columns")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("column classes")
    fig.savefig(path, dpi=150)
    return path


def sites_in_context(annotated: pd.DataFrame, path: Path, distance_label: str) -> Path:
    fig, ax = _figure(size=(7, 4.5))
    for cls in CLASS_ORDER:
        sub = annotated[annotated["class"] == cls]
        if sub.empty:
            continue
        ax.scatter(
            sub["min_distance"], sub["posterior"],
            s=18, alpha=0.75, label=cls, color=PALETTE[cls],
            edgecolors="none",
        )
    ax.set_xlabel(distance_label)
    ax.set_ylabel("posterior Qk")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("divergence against structural context")
    fig.savefig(path, dpi=150)
    return path


def write(outdir: Path, run: RunSummary, classified: pd.DataFrame) -> list[Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = [funnel(run, outdir / "funnel.png"), classes(classified, outdir / "classes.png")]
    (outdir / "summary.txt").write_text(
        run.as_text() + "\n\n" + summary(classified).to_string(index=False) + "\n"
    )
    written.append(outdir / "summary.txt")
    return written
