#!/usr/bin/env bash

usage() { 
    local err=${1:-""};
    cat <<EOF
Usage: $0 [options] 

Required:
    -s FILE   Input snakefile
    -o DIR    Output directory

Optional:
    -c FILE   Snakemake configuration file
    -n        Run snakemake in dry run mode
    -t DIR    The temporary directory to use for Java.
EOF
# shellcheck disable=SC2188
>&2;
    echo -e "\n$err" >&2;
    exit 1;
}

dry_run=""
tmp_dir=""

while getopts "s:o:c:nt:" flag; do
    case "${flag}" in
        s) snakefile=${OPTARG};;
        o) out_dir=${OPTARG};;
        c) config_file=${OPTARG};;
        n) dry_run="-n";;
        t) tmp_dir=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

extra_args=""
if [ -z "${snakefile}" ]; then
    usage "Missing required parameter -s";
fi
if [ -z "${out_dir}" ]; then
    usage "Missing required parameter -o";
fi
if [ -n "${config_file}" ]; then
    extra_args="--configfile $config_file";
fi


# shellcheck disable=SC1091
source "$(dirname "$0")"/common.sh
cores=$(find_core_limit)
mem_gb=$(find_mem_limit_gb)
log "Number of cores: $cores"
log "Memory limit: $mem_gb GB"


# Set the temporary directory for Java executables if necessary.
# These variables must be explicitly set before we run snakemake
if [ "$tmp_dir" != "" ]; then
    export _JAVA_OPTIONS=-Djava.io.tmpdir="$tmp_dir"
    export JAVA_TOOL_OPTIONS=-Djava.io.tmpdir="$tmp_dir"
fi

# Run Snakemake pipeline
set -euo pipefail
# shellcheck disable=SC2086
snakemake \
  --printshellcmds \
  --reason \
  --nocolor \
  --keep-going \
  --rerun-incomplete \
  --jobs "$cores" \
  --resources "mem_gb=$mem_gb" \
  --snakefile $snakefile \
  --directory $out_dir \
  $dry_run \
  $extra_args;

log "All done!"
