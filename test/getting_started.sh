#!/bin/bash
# Test glnexus_cli as suggested by the Getting Started guide on the wiki:
# https://github.com/dnanexus-rnd/GLnexus/wiki/Getting-Started
set -e -o pipefail

echo "--- testing Getting Started examples"

HERE=$(dirname "$0")

if [ ! -f /tmp/dv_platinum6_chr21_gvcf.tar ]; then
    wget -O /tmp/dv_platinum6_chr21_gvcf.tar https://raw.githubusercontent.com/wiki/dnanexus-rnd/GLnexus/data/dv_platinum6_chr21_gvcf.tar
fi

rm -rf /tmp/dv_platinum6_chr21_gvcf/ GLnexus.DB/
tar xf /tmp/dv_platinum6_chr21_gvcf.tar -C /tmp

echo -e "chr21\t0\t48129895" > /tmp/dv_platinum6_chr21_gvcf/hg19_chr21.bed
time $HERE/../glnexus_cli --config DeepVariant --squeeze -a --bed /tmp/dv_platinum6_chr21_gvcf/hg19_chr21.bed \
    /tmp/dv_platinum6_chr21_gvcf/*.gvcf.gz > /tmp/dv_platinum6_chr21_gvcf/dv_platinum6_chr21.bcf

echo "--- Getting Started examples OK"
