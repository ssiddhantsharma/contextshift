"""Gene neighbourhood conservation by group.

FlaGs `_operon.tsv` is positional and unlabelled; columns below are taken from
its writer, where `species` is `<acc>|<species>` and `ids` is `<acc>#<n>`.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from ..partition import Partition
from ..schema import NEIGHBOURS

OPERON_COLUMNS = [
    "species", "length", "query_strand", "neighbour_strand", "cluster",
    "rel_start", "rel_end", "start", "end", "ids",
]
UNCLUSTERED = "0"


def parse_operon_tsv(path: Path) -> pd.DataFrame:
    """Read a FlaGs _operon.tsv into the NEIGHBOURS schema."""
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str).iloc[:, : len(OPERON_COLUMNS)]
    raw.columns = OPERON_COLUMNS
    return _to_neighbours(raw)


def _to_neighbours(raw: pd.DataFrame) -> pd.DataFrame:
    df = raw.copy()
    df["member_id"] = df["species"].str.split("|").str[0]
    df["neighbour_id"] = df["ids"].str.split("#").str[0]
    df["cluster_id"] = df["cluster"].astype("string")
    df["start"] = pd.to_numeric(df["start"], errors="coerce")

    # FlaGs writes no offset column; rank by coordinate with the query at 0
    out = []
    for member, grp in df.sort_values("start").groupby("member_id", sort=False):
        grp = grp.reset_index(drop=True)
        self_rows = grp.index[grp["neighbour_id"] == member].tolist()
        origin = self_rows[0] if self_rows else len(grp) // 2
        grp["offset"] = grp.index - origin
        out.append(grp)

    joined = pd.concat(out, ignore_index=True)
    cols = joined[["member_id", "offset", "neighbour_id", "cluster_id"]].copy()
    cols["offset"] = cols["offset"].astype("Int64")
    return NEIGHBOURS.validate(cols)


def parse_outdesc(path: Path) -> dict[str, str]:
    """Map neighbour accession to product description."""
    out: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            out[parts[1].strip()] = parts[2].strip()
    return out


def annotate(neighbours: pd.DataFrame, descriptions: dict[str, str]) -> pd.DataFrame:
    df = neighbours.copy()
    df["annotation"] = df["neighbour_id"].map(descriptions).astype("string")
    return NEIGHBOURS.validate(df)


def conservation(
    neighbours: pd.DataFrame, partition: Partition, min_fraction: float = 0.5
) -> pd.DataFrame:
    """Fraction of each group's members carrying each neighbour cluster."""
    df = neighbours.copy()
    df["group"] = df["member_id"].map(partition.labels)
    df = df[df["group"].notna()]

    sizes = df.groupby("group")["member_id"].nunique()
    counts = (
        df.groupby(["group", "cluster_id"])["member_id"].nunique().reset_index(name="n_members")
    )
    counts["group_size"] = counts["group"].map(sizes)
    counts["fraction"] = counts["n_members"] / counts["group_size"]
    counts["conserved"] = (counts["fraction"] >= min_fraction) & (
        counts["cluster_id"] != UNCLUSTERED
    )
    return counts.sort_values(["group", "fraction"], ascending=[True, False]).reset_index(drop=True)


def shared_and_private(conservation_table: pd.DataFrame) -> dict[str, object]:
    """Clusters conserved in every group versus in exactly one."""
    cons = conservation_table[conservation_table["conserved"]]
    groups = sorted(conservation_table["group"].unique())
    by_group = {g: set(cons[cons["group"] == g]["cluster_id"]) for g in groups}

    shared = set.intersection(*by_group.values()) if by_group else set()
    private = {
        g: sorted(s - set.union(*(by_group[o] for o in groups if o != g), set()))
        for g, s in by_group.items()
    }
    return {
        "groups": groups,
        "n_conserved": {g: len(s) for g, s in by_group.items()},
        "shared": sorted(shared),
        "private": private,
        "groups_with_none": [g for g, s in by_group.items() if not s],
    }


# --- gene neighbourhoods from a GFF ----------------------------------------

GFF_ID = re.compile(r"(?:^|;)(?:protein_id|ID)=(?:cds-)?([^;]+)")
GFF_PRODUCT = re.compile(r"(?:^|;)product=([^;]+)")


def read_gff_cds(path: Path) -> pd.DataFrame:
    """Coding features from a GFF, ordered along each sequence."""
    rows = []
    for line in Path(path).read_text(errors="replace").splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        ident = GFF_ID.search(f[8])
        product = GFF_PRODUCT.search(f[8])
        rows.append({
            "contig": f[0],
            "start": int(f[3]),
            "end": int(f[4]),
            "strand": f[6],
            "member_id": ident.group(1) if ident else "",
            "product": product.group(1).replace("%2C", ",") if product else "",
        })
    df = pd.DataFrame(rows).drop_duplicates(subset=["contig", "start", "end", "member_id"])
    df = df.sort_values(["contig", "start"]).reset_index(drop=True)
    df["index_on_contig"] = df.groupby("contig").cumcount()
    return df


def neighbourhood(
    cds: pd.DataFrame, member_id: str, window: int = 5
) -> tuple[pd.DataFrame, bool]:
    """Genes within `window` positions of one member, and whether the window
    ran off the end of the sequence.

    The flag matters: a gene near a contig edge has an incomplete
    neighbourhood, so concluding anything from the absence of a neighbour
    there is unsound.
    """
    hit = cds[cds["member_id"] == member_id]
    if hit.empty:
        raise KeyError(member_id)
    row = hit.iloc[0]
    same = cds[cds["contig"] == row["contig"]]
    centre = int(row["index_on_contig"])
    lo, hi = centre - window, centre + window
    truncated = lo < 0 or hi > int(same["index_on_contig"].max())
    near = same[(same["index_on_contig"] >= lo) & (same["index_on_contig"] <= hi)].copy()
    near["offset"] = near["index_on_contig"] - centre
    return near.reset_index(drop=True), truncated


def classify_by_neighbours(
    cds: pd.DataFrame,
    members: list[str],
    pattern: str,
    window: int = 5,
) -> pd.DataFrame:
    """Label each member by whether a neighbour's product matches `pattern`.

    `edge` marks members whose window was truncated by a sequence end; a
    negative call there is not safe to trust.
    """
    rx = re.compile(pattern, re.I)
    rows = []
    for member in members:
        try:
            near, truncated = neighbourhood(cds, member, window)
        except KeyError:
            rows.append({"member_id": member, "label": "not_in_gff", "n_matching": 0,
                         "edge": False, "nearest_offset": pd.NA, "matches": ""})
            continue
        others = near[near["member_id"] != member]
        # astype(bool) matters: an empty .map() result is object-dtype, which
        # pandas reads as column selection rather than a mask
        mask = others["product"].map(lambda p: bool(rx.search(p))).astype(bool)
        matched = others[mask]
        nearest = (matched["offset"].abs().min() if len(matched) else pd.NA)
        rows.append({
            "member_id": member,
            "label": "associated" if len(matched) else "solo",
            "n_matching": len(matched),
            "edge": truncated,
            "nearest_offset": nearest,
            "matches": "; ".join(sorted(set(matched["product"]))[:3]),
        })
    return pd.DataFrame(rows)
