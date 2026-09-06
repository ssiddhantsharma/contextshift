"""Report figures render and never hide an exclusion."""

import pandas as pd
import pytest

from contextshift import report
from contextshift.join import CORE, DETERMINANT, classify
from contextshift.partition import Partition, Provenance


def classified_frame():
    sites, cons = [], []
    for col, post in [(10, 0.05), (20, 0.97)]:
        sites.append({
            "family": "F", "partition": "context", "group_a": "A", "group_b": "B",
            "column": col, "test": "type2", "posterior": post,
        })
        rate_all = 0.05 if post < 0.5 else 1.2
        for scope, rate in [("all", rate_all), ("A", 0.1), ("B", 0.1)]:
            cons.append({"family": "F", "column": col, "scope": scope, "rate": rate})
    return classify(pd.DataFrame(sites), pd.DataFrame(cons))


def test_summary_lists_dropped_counts():
    run = report.RunSummary()
    run.stage("loci", 100, dropped=5)
    run.stage("members", 300)
    text = run.as_text()
    assert "-5 dropped" in text
    assert "300" in text


def test_notes_appear_in_the_summary():
    run = report.RunSummary(notes=["type-I ran through a shim"])
    run.stage("sites", 10)
    assert "shim" in run.as_text()


def figures(tmp_path):
    run = report.RunSummary()
    run.stage("loci", 100, dropped=5)
    run.stage("members", 280, dropped=20)
    partition = Partition(
        name="subtype",
        labels={f"m{i}": ("clean" if i < 30 else "messy") for i in range(60)},
        provenance=Provenance(source="t"),
        caveats={"messy": "not a clade"},
    )
    df = classified_frame()
    df["posterior"] = [0.05, 0.97]
    df["min_distance"] = [12.0, 3.4]
    return {
        "funnel": lambda p: report.funnel(run, p),
        "power": lambda p: report.power(partition, p),
        "classes": lambda p: report.classes(classified_frame(), p),
        "sites": lambda p: report.sites_in_context(df, p, "distance to DNA (A)"),
    }


@pytest.mark.parametrize("name", ["funnel", "power", "classes", "sites"])
def test_panel_renders(tmp_path, name):
    out = figures(tmp_path)[name](tmp_path / f"{name}.png")
    assert out.exists() and out.stat().st_size > 1000


def test_write_emits_figures_and_a_summary(tmp_path):
    run = report.RunSummary()
    run.stage("columns", 2)
    written = report.write(tmp_path, run, classified_frame())
    assert {p.name for p in written} == {"funnel.png", "classes.png", "summary.txt"}
    text = (tmp_path / "summary.txt").read_text()
    assert CORE in text and DETERMINANT in text
