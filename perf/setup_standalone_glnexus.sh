#!/bin/bash

set -ex -o pipefail

source ~/dx-toolkit/environment

# debug packages
sudo apt-get install libjemalloc1-dbg

# download data sets
echo "Downloading data sets"
mkdir -p /mnt/U/data
cd /mnt/U/data
if [ ! -e "vcr_b37.bed" ]; then
    dx download RGC_Coriell_NODOWNLOAD:/
    cp vcr_b37.bed ranges.bed
fi
if [ ! -e "gvcfs1024.tar" ]; then
    mkdir gvcf
    pushd gvcf
    dx download RGC_Coriell_NODOWNLOAD:/Fake/gvcfs1024.tar
    tar xvf gvcfs1024.tar --strip-components=1
    popd
fi

find gvcf -type f > all_gvcfs.txt
wc -l all_gvcfs.txt

echo "Initilize database"
glnexus_cli init GLnexus.db $(find gvcf -type f | head -n 1)

echo "Loading data"
cat all_gvcfs.txt | time glnexus_cli load --and-delete GLnexus.db -
