"""Shelling out to external tools, with the version recorded."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


class ToolMissing(RuntimeError):
    pass


@dataclass(frozen=True)
class Tool:
    name: str
    version_args: tuple[str, ...] = ("--version",)
    install_hint: str = ""

    @property
    def path(self) -> str | None:
        return shutil.which(self.name)

    def require(self) -> str:
        p = self.path
        if p is None:
            hint = f" ({self.install_hint})" if self.install_hint else ""
            raise ToolMissing(f"{self.name} not found on PATH{hint}")
        return p

    def version(self) -> str:
        if self.path is None:
            return "missing"
        try:
            out = subprocess.run(
                [self.name, *self.version_args],
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
            return (out.stdout or out.stderr).strip().splitlines()[0]
        except Exception:
            return "unknown"


MMSEQS = Tool("mmseqs", ("version",), "conda install -c bioconda mmseqs2")
MAFFT = Tool("mafft", ("--version",), "conda install -c bioconda mafft")
FOLDMASON = Tool("foldmason", ("version",), "conda install -c bioconda foldmason")
IQTREE = Tool("iqtree2", ("--version",), "conda install -c bioconda iqtree")
RATE4SITE = Tool("rate4site", ("-h",), "conda install -c bioconda rate4site")
MEME = Tool("meme", ("-version",), "conda install -c bioconda meme")
CCTYPER = Tool("cctyper", ("--version",), "pip install cctyper")


def run(cmd: list[str], cwd: Path | None = None, timeout: int = 3600) -> subprocess.CompletedProcess:
    result = subprocess.run(
        cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout, check=False
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"{cmd[0]} exited {result.returncode}\ncmd: {' '.join(cmd)}\n{result.stderr[-2000:]}"
        )
    return result


def toolchain_report() -> str:
    tools = [MMSEQS, MAFFT, FOLDMASON, IQTREE, RATE4SITE, MEME, CCTYPER]
    width = max(len(t.name) for t in tools)
    lines = []
    for t in tools:
        mark = "x" if t.path else " "
        lines.append(f"  [{mark}] {t.name:{width}s}  {t.version()}")
    return "\n".join(lines)
