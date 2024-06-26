## base container
FROM mambaorg/micromamba:0.15.3 as base_container
USER root
RUN apt-get update \
    && apt-get install --yes rename tini procps curl \
    && apt-get autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

USER micromamba
ENTRYPOINT ["tini", "--"]
CMD ["/bin/bash"]

## main container
ARG CONDA_FILE=containers/environment.yml
FROM base_container
# adding opencontainer labels to link registry to github repository
LABEL org.opencontainers.image.title="openmm-nf"
LABEL org.opencontainers.image.description="A simple workflow to evaluate protein stability with respect to the wildtype"
LABEL org.opencontainers.image.url="https://github.com/stracquadaniolab/openmm-nf"
LABEL org.opencontainers.image.documentation="https://github.com/stracquadaniolab/openmm-nf"
LABEL org.opencontainers.image.source="https://github.com/stracquadaniolab/openmm-nf"
LABEL org.opencontainers.image.vendor="stracquadaniolab"
LABEL org.opencontainers.image.authors="Josh David Littlefair"
LABEL org.opencontainers.image.revision="v0.1.57"

COPY --chown=micromamba:micromamba ${CONDA_FILE} /tmp
RUN micromamba install -y -n base -f /tmp/`basename ${CONDA_FILE}` && \
    micromamba clean --all --yes
