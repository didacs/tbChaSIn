FROM mambaorg/micromamba:1.5.6 as base
# Install system-level dependencies as root
USER root
RUN apt-get update && apt-get install -y procps
USER $MAMBA_USER
# Set working directory
WORKDIR /tbChaSIn
# Lock files can cause problems
#RUN micromamba config set use_lockfiles False
# Downloads can be slow sometimes
#RUN micromamba config set remote_connect_timeout_secs 300
#RUN micromamba config set remote_backoff_factor 10
#RUN micromamba config set remote_max_retries 5
# This completely disables timeouts due to slow downloads but it can result in a hung install
#ENV MAMBA_NO_LOW_SPEED_LIMIT=1

FROM base as fastp
COPY --chown=$MAMBA_USER:$MAMBA_USER ./mamba/fastp.yml .
RUN micromamba install -y -n base -f fastp.yml \
    && micromamba clean --all --yes \
    && rm -f fastp.yml
ARG MAMBA_DOCKERFILE_ACTIVATE=1

FROM base as fastqc
COPY --chown=$MAMBA_USER:$MAMBA_USER ./mamba/fastqc.yml .
RUN micromamba install -y -n base -f fastqc.yml \
    && micromamba clean --all --yes \
    && rm -f fastqc.yml
ARG MAMBA_DOCKERFILE_ACTIVATE=1

FROM base as multiqc
COPY --chown=$MAMBA_USER:$MAMBA_USER ./mamba/multiqc.yml .
RUN micromamba install -y -n base -f multiqc.yml \
    && micromamba clean --all --yes \
    && rm -f multiqc.yml
ARG MAMBA_DOCKERFILE_ACTIVATE=1

FROM base as picard
COPY --chown=$MAMBA_USER:$MAMBA_USER ./mamba/picard.yml .
RUN micromamba install -vvv -y -n base -f picard.yml \
    && micromamba clean --all --yes \
    && rm -f picard.yml
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# TODO: processes actually use three different combinations of tools. Ideally, we
# would have a separate environment.yml file for each and compose them using
# `micromamba env update -n base -f <tool>.yml`, but it doesn't look like env update
# will be available until micromamba 2.0.
# 1. fgbio
# 2. fgbio + samtools
# 3. fgbio + samtools + bwa
FROM base as fgbio
COPY --chown=$MAMBA_USER:$MAMBA_USER ./mamba/fgbio.yml .
RUN micromamba install -y -n base -f fgbio.yml \
    && micromamba clean --all --yes \
    && rm -f fgbio.yml
ARG MAMBA_DOCKERFILE_ACTIVATE=1

FROM base as tools-base
# Download the Aurora RDS root certificate bundle to the default location. This is required to
# connect to the Benchling Warehouse database.
# See https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/UsingWithRDS.SSL.html
ADD --chown=$MAMBA_USER:$MAMBA_USER \
    https://truststore.pki.rds.amazonaws.com/global/global-bundle.pem \
    /home/$MAMBA_USER/.postgresql/root.crt
# Set environment variables to require root CA verification by default and make these settings 
# apply to new shells. Use .zshrc for ZSH.
ENV PGSSLMODE="verify-ca"

# Build tomebio-tools
FROM tools-base as tools-build
USER root
RUN apt-get install -y build-essential
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER ./mamba/dev.yml .
RUN micromamba install -y -n base -f dev.yml \
    && micromamba clean --all --yes \
    && rm -Rf dev.yml 
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ADD --chown=$MAMBA_USER:$MAMBA_USER https://install.python-poetry.org install-poetry.py
RUN python3 install-poetry.py
ENV PATH="${PATH}:~/.local/bin"
# Install dependencies
COPY --chown=$MAMBA_USER:$MAMBA_USER ./src/python/pyproject.toml .
COPY --chown=$MAMBA_USER:$MAMBA_USER ./src/python/poetry.lock .
RUN poetry install --no-root
# pybedtools does not provide a wheel, so we compile it here to avoid having
# gcc installed in the final image
RUN export pybedtools_version=$(poetry show pybedtools \
    | grep version \
    | sed -rn 's/.*: (.+)/\1/p' \
    | tr -d '[:space:]' ) \
    && pip download --no-binary=:all: --no-deps pybedtools=="${pybedtools_version}" \
    && pip wheel --no-deps -w dist "pybedtools-${pybedtools_version}.tar.gz"
# Build the wheel - do this last so we only have to re-run this
# part when the code changes 
COPY --chown=$MAMBA_USER:$MAMBA_USER ./src/python/ .
RUN poetry build

FROM tools-base as tools
COPY --chown=$MAMBA_USER:$MAMBA_USER ./mamba/tools.yml .
RUN micromamba install -y -n base -f tools.yml \
    && micromamba clean --all --yes \
    && rm tools.yml
ARG MAMBA_DOCKERFILE_ACTIVATE=1
# Install pytomebio
COPY --from=tools-build --chown=$MAMBA_USER:$MAMBA_USER /tbChaSIn/dist/*.whl .
RUN pip install *.whl \
    && rm -Rf *.whl
