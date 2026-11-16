#!env bash
set -o pipefail
wget -c -t 5 -w 5 -T 60 http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.2.tar.gz
tar -xzvf eukcc2_db_ver_1.2.tar.gz -C $CONDA_PREFIX
rm eukcc2_db_ver_1.2.tar.gz
chmod -R 777 $CONDA_PREFIX/eukcc2_db_ver_1.2
