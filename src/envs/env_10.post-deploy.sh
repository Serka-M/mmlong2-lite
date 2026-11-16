#!env bash
set -o pipefail
git clone https://github.com/LorMeBioAI/LorBin $CONDA_PREFIX/LorBin
pip3 install --no-cache-dir $CONDA_PREFIX/LorBin/dist/lorbin-*.tar.gz
