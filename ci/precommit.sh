#!/usr/bin/env bash

###############################################################################
# Script that should be run pre-commit after making any changes to the pytomebio
# package / subdirectory.
#
# Runs:
#   Unit tests
#   Linting
#   Type checking
#   Style checking
###############################################################################

set -e

failures=""

usage() {
    cat <<EOF
Usage: $0 [options] 

Optional:
    -f    Automatically fix formatting.
EOF
    # shellcheck disable=SC2188
    >&2
    exit 1
}

function banner() {
    echo
    echo "================================================================================"
    echo "$*"
    echo "================================================================================"
    echo
}

#####################################################################
# Takes two parameters, a "name" and a "command".
# Runs the command and prints out whether it succeeded or failed, and
# also tracks a list of failed steps in $failures.
#####################################################################
function run() {
    local name=$1
    local cmd=$2

    banner "Running $name [$cmd]"
    set +e
    $cmd
    exit_code=$?
    set -e

    if [[ $exit_code == 0 ]]; then
        echo "Passed $name"
    else
        echo "Failed $name [$cmd]"
        if [ -z "$failures" ]; then
            failures="$failures $name"
        else
            failures="$failures, $name"
        fi
    fi
}

fix='false'
while getopts "f" flag; do
    case "${flag}" in
    f) fix='true' ;;
    *) usage ;;
    esac
done
shift $((OPTIND - 1))

parent=$(cd "$(dirname "$0")" && pwd -P)
repo_root="$(dirname "${parent}")"
root=${repo_root}/src/python

if [[ -z ${CONDA_DEFAULT_ENV} ]]; then
    banner "Conda not active. pytomebio conda environment must be active."
    exit 1
fi

if ${fix}; then
    pushd "$root" >/dev/null
    banner "Executing in conda environment ${CONDA_DEFAULT_ENV} in directory ${root}"
    run "Snakemake Style Formatting" "snakefmt --line-length 99 ${repo_root}/src/snakemake"
    run "Python Style Formatting" "ruff format --config=$parent/ruff.toml pytomebio"
    popd >/dev/null
fi

pushd "$root" >/dev/null
banner "Executing in conda environment ${CONDA_DEFAULT_ENV} in directory ${root}"
run "Shell Check (precommit)" "shellcheck ${repo_root}/ci/precommit.sh"
run "Shell Check (src/scripts)" "shellcheck ${repo_root}/src/scripts/*sh"
run "Snakemake Style Checking" "snakefmt --check --line-length 99 ${repo_root}/src/snakemake"
run "Python Style Checking" "ruff format --config=$parent/ruff.toml --check pytomebio"
run "Python Linting" "ruff check --config=$parent/ruff.toml pytomebio"
run "Python Type Checking" "mypy -p pytomebio --config $parent/mypy.ini"
run "Python Unit Tests" "pytest -vv -r sx pytomebio"
popd >/dev/null

if [ -z "$failures" ]; then
    banner "Precommit Passed"
else
    banner "Precommit Failed with failures in: $failures"
    exit 1
fi
