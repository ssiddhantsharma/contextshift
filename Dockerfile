# contextshift with its optional tools, including a working DIVERGE build.
#
# DIVERGE cannot be pip-installed: its sdist omits a requirements.txt that
# setup.py reads, and its default src/ tree is MSVC-only. This builds from
# src_linux and drops a `version` macro that collides with a pybind11 member.
# Upstream fixes: zjupgx/diverge4 PRs #8, #9, #10.

FROM mambaorg/micromamba:1.5-jammy

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 mmseqs2 mafft iqtree meme \
    && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /opt
RUN git clone --depth 1 https://github.com/zjupgx/diverge4.git \
    && cd diverge4 \
    && python - <<'PY'
import re
from pathlib import Path

# build the portable tree, not the MSVC one
s = Path("setup.py").read_text()
s = s.replace("glob(f'src/{module}/*.c*')", "glob(f'src_linux/{module}/*.c*')")
Path("setup.py").write_text(s)

# `#define version "V1.0"` collides with a member named version in
# pybind11 detail/internals.h; the macro is unused
target = b'#define version "V1.0"\n'
for p in Path("src_linux").rglob("common.h"):
    b = p.read_bytes()
    if target in b:
        p.write_bytes(b.replace(target, b""))

# Gu99/Rvs/TypeOneAnalysis call an undefined get_colnames (upstream PR #8)
b = Path("diverge/binding.py").read_text()
if "def get_colnames" not in b:
    b = b.replace("def load_tree_file(",
                  "def get_colnames(names):\n    return list(names)\n\n\ndef load_tree_file(", 1)
    Path("diverge/binding.py").write_text(b)
PY
RUN pip install --no-cache-dir pybind11 setuptools wheel \
    && pip install --no-cache-dir --no-build-isolation ./diverge4 \
    && rm -rf /opt/diverge4

WORKDIR /work
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/contextshift
RUN pip install --no-cache-dir /opt/contextshift

ENTRYPOINT ["contextshift"]
CMD ["doctor"]
