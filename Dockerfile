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

# CCTyper is absent on arm64. Not because a dependency lacks an ARM build --
# every one of them has one -- but because the bioconda recipe pins
# `prodigal >=2.0,<=2.6.2` and only prodigal 2.6.3 ships linux-aarch64.
# Add `cctyper` to the list below on x86_64, or once bioconda-recipes#68871 lands.
# The `loci` stage only parses typer output, so typing can be done elsewhere.
USER $MAMBA_USER
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 mmseqs2 mafft iqtree meme \
    && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root
RUN mkdir -p /build /work && chown -R $MAMBA_USER /build /work
USER $MAMBA_USER

WORKDIR /build
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/patch_diverge.py .
RUN git clone --depth 1 https://github.com/zjupgx/diverge4.git \
    && python patch_diverge.py diverge4 \
    && pip install --no-cache-dir pybind11 setuptools wheel \
    && pip install --no-cache-dir --no-build-isolation ./diverge4 \
    && rm -rf /build/diverge4

WORKDIR /work
COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp/contextshift
RUN pip install --no-cache-dir /tmp/contextshift && rm -rf /tmp/contextshift

# MAMBA_DOCKERFILE_ACTIVATE only affects build-time RUN steps; the base image's
# entrypoint activates the environment for the running container.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "contextshift"]
CMD ["doctor"]
