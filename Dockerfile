# This Dockerfile copies the entire codebase into an image based on 
# micromamba and creates the environment with the name 'base' (as recommended
# here: https://micromamba-docker.readthedocs.io/en/latest/quick_start.html).
# To build this image run:
# docker build [--platform linux/amd64] -t tomebio/cryptic-seq:<tag>
# the '--platform' argument is only required on a Mac with Apple silicon.
# If there are issues and you need to debug:
# 1. Add the verbose option (-vvv) to micromamba install
# 2. Create a docker build environment with a larger log size:
# docker buildx create --use --name larger_log --driver-opt env.BUILDKIT_STEP_LOG_MAX_SIZE=50000000
# 3. Build using the environment:
# docker buildx build [--platform linux/amd64] --progress plain --no-cache --load -t tomebio/cryptic-seq:<tag> .
FROM mambaorg/micromamba:1.5.6 as base
# Install system-level dependencies as root
USER root
RUN apt-get update \
    && apt-get install -y procps
USER $MAMBA_USER
# Set working directory to where we will copy the repo 
WORKDIR /tbChaSIn
# Copy environment into the container
COPY --chown=$MAMBA_USER:$MAMBA_USER ./environment.yml .
# Lock files can cause problems
RUN micromamba config set use_lockfiles False
# Downloads can be slow sometimes
#RUN micromamba config set remote_connect_timeout_secs 300
#RUN micromamba config set remote_backoff_factor 10
#RUN micromamba config set remote_max_retries 5
# This completely disables timeouts due to slow downloads but it can result in a hung install
#ENV MAMBA_NO_LOW_SPEED_LIMIT=1
# Install packages from conda environment (*must* use name 'base')
RUN micromamba install -y -n base -f environment.yml
# Clean up
RUN micromamba clean --all --yes
# Activate the environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1
# Copy entire project into the container
# We do this in a separate step so we don't have to re-install the
# environment any time there is a change to the code.
# TODO: add '--exclude=environment.yml' once that option becomes stable
COPY --chown=$MAMBA_USER:$MAMBA_USER . .

# Build pytomebio wheel
FROM base as build
# Install poetry
ADD --chown=$MAMBA_USER:$MAMBA_USER https://install.python-poetry.org install-poetry.py
RUN python3 install-poetry.py
ENV PATH="${PATH}:~/.local/bin"
# Build pytomebio
RUN poetry build

# Create file image by installing wheel and creating pipeline wrappers
FROM base
COPY --from=build --chown=$MAMBA_USER:$MAMBA_USER /tbChaSIn/dist/*.whl ./dist/
RUN pip install dist/*.whl \
    && rm -Rf dist
# Download the Aurora RDS root certificate bundle to the default location. This is required to
# connect to the Benchling Warehouse database.
# See https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/UsingWithRDS.SSL.html
ADD --chown=$MAMBA_USER:$MAMBA_USER \
    https://truststore.pki.rds.amazonaws.com/global/global-bundle.pem \
    /home/$MAMBA_USER/.postgresql/root.crt
# Set environment variables to require root CA verification by default and make these settings 
# apply to new shells. Use .zshrc for ZSH.
ENV PGSSLMODE="verify-ca"
# Run tests
RUN bash ci/precommit.sh
# Wrap pipeline(s)
RUN snk install src/snakemake/cryptic-seq
