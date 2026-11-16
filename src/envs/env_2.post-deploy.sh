#!env bash
set -o pipefail
dir=$(pwd)
cd $CONDA_PREFIX
git clone https://github.com/bluenote-1577/rust-anvio-mis.git rust-anvio-mis
cd rust-anvio-mis
cargo generate-lockfile
cargo update -p hts-sys --precise 2.1.1
cargo install --locked --path . --root $CONDA_PREFIX
cd ..
yes | rm -r rust-anvio-mis
cd $dir
