FROM mambaorg/micromamba:1.5.8

WORKDIR /ara

COPY --chown=$MAMBA_USER:$MAMBA_USER ./requirements.yaml /tmp/

COPY --chown=$MAMBA_USER:$MAMBA_USER . .

RUN micromamba install -y -n base -f /tmp/requirements.yaml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN perl setup.pl

ENV PATH="/opt/conda/bin:/ara:$PATH"

CMD perl ara.pl -h