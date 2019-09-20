#!/bin/bash

# prepare the dataset staged at
#   https://raw.githubusercontent.com/wiki/dnanexus-rnd/GLnexus/data/dv_1000G_ALDH2_gvcf.tar
# containing  the DeepVariant gVCF slices for ALDH2 locus from the 1000 Genomes Project data.
# This is used in the GLnexus "Getting Started" tutorial and integration tests.
#
# dependencies:
# - gsutil
# - parallel
# - tabix (& bgzip)
# - ifne (apt-get install moreutls)

set -Eeuo pipefail
rm -rf dv_1000G_ALDH2_gvcf/ dv_1000G_ALDH2_gvcf.tar
mkdir dv_1000G_ALDH2_gvcf
cd dv_1000G_ALDH2_gvcf

gsutil ls gs://brain-genomics-public/research/cohort/1000G/dv_vcf/*.g.vcf.gz \
    | parallel --retries 3 -t 'set -Eeuo pipefail; tabix -h {} chr12:111,760,000-111,820,000 | ifne -n sh -c "exit 99" | bgzip -c > $(basename {})'
# ifne used because tabix with gs:// is prone to producing empty output with 0 exit status in case
# of sporadic errors

cd ..
tar -cvf dv_1000G_ALDH2_gvcf.tar dv_1000G_ALDH2_gvcf/*.g.vcf.gz
