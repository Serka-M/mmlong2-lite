#!env bash
set -o pipefail
wget https://sourceforge.net/projects/bbmap/files/BBMap_40.02.tar.gz/download -O BBTools.tar.gz
tar -xzf BBTools.tar.gz -C $CONDA_PREFIX
rm BBTools.tar.gz
ln -s $CONDA_PREFIX/bbmap/quickbin.sh $CONDA_PREFIX/bin/quickbin.sh
chmod +x $CONDA_PREFIX/bin/quickbin.sh
