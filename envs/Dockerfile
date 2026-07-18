# Default (GenMap) environment for reproducible builds and Zenodo provenance.
FROM mambaorg/micromamba:1.5-jammy
COPY environment.genmap.yml /tmp/environment.genmap.yml
RUN micromamba install -y -n base -f /tmp/environment.genmap.yml && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
WORKDIR /work
