#!/usr/bin/env bash

###############################################################################
# Script that builds Docker images. For nextflow, it is expected that there is
# a 1:1 correspondence between Mamba config files in mamba/ and Docker
# images. Images are tagged with 'tomebio/<conf>:<version>', where <conf> is
# the name of the Mamba config file without the .yml extension, and <version>
# defaults to 1.0. For snakemake, a single image is built with all the
# dependencies, and it is tagged with 'tomebio/cryptic-seq:<version>'.
###############################################################################

set -e

usage() {
    cat <<EOF
Usage: $0 [options] [name...]

All images are built unless specific images are listed by name.
Defaults to tagging images with "dev" ECR repo uri

Optional:
    -s          Build the 'all-in-one' image for the snakemake pipeline
    -v VERSION  Version to use when tagging images
    -f          Force Docker to rebuild all images
    -m          Build image to run on ARM Mac
    -p          Push images to ECR
    -x          Tag images with "prod" ECR repo URI
EOF
    # shellcheck disable=SC2188
    >&2
    exit 1
}

parent=$(cd "$(dirname "$0")" && pwd -P)
root="$(dirname "${parent}")"
confs="${root}/mamba"

# Source variables from pipeline.env
# shellcheck disable=SC1091
source "${root}/src/nextflow/cryptic-seq/pipeline.env"
# shellcheck disable=SC2153
version="$PIPELINE_VERSION"
# shellcheck disable=SC2153
aws_region="$AWS_REGION"

snakemake=false
rosetta_options=()
force_option=""
push=false
prod=false
while getopts "v:smfpx" flag; do
    case "${flag}" in
    s) snakemake=true ;;
    v) version=${OPTARG} ;;
    m) rosetta_options=(--platform linux/amd64 --load) ;;
    f) force_option="--no-cache" ;;
    p) push=true ;;
    x) prod=true ;;
    *) usage ;;
    esac
done
shift $((OPTIND - 1))

if ${prod}; then
    ecr_repo_root="$PROD_ECR_REPO_ROOT"
else
    ecr_repo_root="$DEV_ECR_REPO_ROOT"
fi

build_target() {
    docker build \
        --progress=plain \
        "${rosetta_options[@]}" \
        ${force_option} \
        --target "${1}" \
        -t "${ecr_repo_root}/tbchasin-${1}:${version}" \
        -f "${root}/docker/Dockerfile" \
        "${root}"
}

build_snakemake() {
    docker build \
        --progress=plain \
        "${rosetta_options[@]}" \
        ${force_option} \
        -t "${ecr_repo_root}/cryptic-seq:${version}" \
        -f "${root}/docker/Dockerfile.snakemake" \
        "${root}"
}

push_target() {
    docker image push "${ecr_repo_root}/tbchasin-${1}:${version}"
}

push_snakemake() {
    docker image push "${ecr_repo_root}/tbchasin-snakemake:${version}"
}

docker_login() {
    aws ecr get-login-password --region "$aws_region" | docker login --username AWS --password-stdin "${ecr_repo_root}"
}

if ${push}; then
    docker_login
fi

if ${snakemake}; then
    build_snakemake
    if ${push}; then
        push_snakemake
    fi
elif [ "$#" -eq 0 ]; then
    # build all images
    for conf_yml in "${confs}"/*.yml; do
        filename=$(basename -- "${conf_yml}")
        target="${filename%.yml}"
        # ignore dev.yml
        if [ "${target}" = "dev" ]; then
            continue
        fi
        build_target "${target}"
        if ${push}; then
            push_target "${target}"
        fi
    done
else
    for target in "$@"; do
        build_target "${target}"
        if ${push}; then
            push_target "${target}"
        fi
    done
fi
