"""Published claims expressed as executable checks.

A claim is a statement taken from a paper, restated so the pipeline can agree
or disagree with it. Disagreement is a result, not an error, so a failing claim
never raises: it is recorded with its observed value.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class Claim:
    name: str
    statement: str
    source: str
    check: Callable[[Any], tuple[bool, float | str]]
    tolerance: str = ""


@dataclass(frozen=True)
class ClaimResult:
    claim: Claim
    upheld: bool
    observed: float | str
    error: str = ""

    @property
    def verdict(self) -> str:
        if self.error:
            return "ERROR"
        return "UPHELD" if self.upheld else "REFUTED"


def evaluate(claims: list[Claim], context: Any) -> list[ClaimResult]:
    results: list[ClaimResult] = []
    for claim in claims:
        try:
            upheld, observed = claim.check(context)
            results.append(ClaimResult(claim, upheld, observed))
        except Exception as exc:  # a claim that cannot be tested is not a claim that passed
            results.append(ClaimResult(claim, False, "", error=str(exc)))
    return results


def report(results: list[ClaimResult]) -> str:
    lines = []
    for r in results:
        lines.append(f"[{r.verdict:7s}] {r.claim.name}")
        lines.append(f"           {r.claim.statement}")
        lines.append(f"           source: {r.claim.source}")
        if r.error:
            lines.append(f"           error: {r.error}")
        else:
            lines.append(f"           observed: {r.observed}")
    return "\n".join(lines)
