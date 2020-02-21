#!/bin/bash
# Test glnexus_cli as suggested by the Getting Started guide on the wiki:
# https://github.com/dnanexus-rnd/GLnexus/wiki/Getting-Started
set -e -o pipefail

echo "--- testing Getting Started examples"

HERE=$(dirname "$0")

if [ ! -f /tmp/dv_1000G_ALDH2_gvcf.tar ]; then
    wget -O /tmp/dv_1000G_ALDH2_gvcf.tar https://raw.githubusercontent.com/wiki/dnanexus-rnd/GLnexus/data/dv_1000G_ALDH2_gvcf.tar
fi

rm -rf /tmp/dv_1000G_ALDH2_gvcf/ GLnexus.DB/
tar xf /tmp/dv_1000G_ALDH2_gvcf.tar -C /tmp

echo -e "chr12\t111760000\t111820000" > /tmp/dv_1000G_ALDH2_gvcf/ALDH2.bed
time $HERE/../glnexus_cli --config DeepVariant --squeeze --bed /tmp/dv_1000G_ALDH2_gvcf/ALDH2.bed \
    /tmp/dv_1000G_ALDH2_gvcf/*.g.vcf.gz > /tmp/dv_1000G_ALDH2_gvcf/dv_1000G_ALDH2.bcf

echo "--- Getting Started examples OK"
