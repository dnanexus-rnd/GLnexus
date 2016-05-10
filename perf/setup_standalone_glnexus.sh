#!/bin/bash

set -ex -o pipefail

source ~/dx-toolkit/environment

# download data sets
echo "Downloading data sets"
mkdir -p /mnt/U/data
cd /mnt/U/data
dx download RGC_Coriell_NODOWNLOAD:/vcr_b37.bed
dx download RGC_Coriell_NODOWNLOAD:/Fake/

tar xvf gvcfs1024.tar --strip-components=1
find gvcf -type f > all_gvcfs.txt
wc -l all_gvcfs.txt
cp vcr_b37.bed ranges.bed

echo "Initilize database"
glnexus_cli init GLnexus.db $(find gvcf -type f | head -n 1)

echo "Loading data"
cat all_gvcfs.txt | time glnexus_cli load --and-delete GLnexus.db -