#!env bash
set -o pipefail
checkm2 database --download --path $CONDA_PREFIX || true
