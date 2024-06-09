FROM mambaorg/micromamba:1.5.8

WORKDIR /ara

COPY --chown=$MAMBA_USER:$MAMBA_USER . .

RUN micromamba install -y -n base -f /ara/requirements.yaml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER ./container/user-settings.mkfg /home/mambauser/.ncbi/user-settings.mkfg

RUN vdb-config --report-cloud-identity yes

RUN perl setup.pl

ENV PATH="/opt/conda/bin:/ara:$PATH"

CMD perl ara.pl -h