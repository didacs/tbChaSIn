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

Optional:
    -s          Build the 'all-in-one' image for the snakemake pipeline
    -v VERSION  Version to use when tagging images
    -f          Force Docker to rebuild all images
    -m          Build image to run on ARM Mac
EOF
    # shellcheck disable=SC2188
    >&2
    exit 1
}

parent=$(cd "$(dirname "$0")" && pwd -P)
root="$(dirname "${parent}")"

project="tomebio"
version="1.0"
snakemake=false
rosetta_options=()
force_option=""
while getopts "v:smf" flag; do
    case "${flag}" in
    s) snakemake=true ;;
    v) version=${OPTARG} ;;
    m) rosetta_options=(--platform linux/amd64 --load) ;;
    f) force_option="--no-cache" ;;
    *) usage ;;
    esac
done
shift $((OPTIND - 1))

build_target() {
    docker build \
        --progress=plain \
        "${rosetta_options[@]}" \
        ${force_option} \
        --target "${1}" \
        -t "${project}/${1}:${version}" \
        -f "${root}/docker/Dockerfile" \
        "${root}"
}

build_snakemake() {
    docker build \
        --progress=plain \
        "${rosetta_options[@]}" \
        ${force_option} \
        -t "${project}/cryptic-seq:${version}" \
        -f "${root}/docker/Dockerfile.snakemake" \
        "${root}"
}

if ${snakemake}; then
    build_snakemake
elif [ "$#" -eq 0 ]; then
    # build all images
    while read -r target; do
        build_target "${target}"
    done <"${root}/docker/.targets"
else
    for target in "$@"; do
        build_target "${target}"
    done
fi
