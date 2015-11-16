#!/bin/bash
# Launch medium-sized tests on DNAnexus. Depends on secure environment
# variables DX_AUTH_TOKEN and GLNEXUS_TEST_PROJECT_ID

set -e -o pipefail

# install dx-toolkit
wget https://wiki.dnanexus.com/images/files/dx-toolkit-current-ubuntu-12.04-amd64.tar.gz
tar zxf dx-toolkit-current-ubuntu-12.04-amd64.tar.gz
source dx-toolkit/environment

# build GLnexus applet
dx cd ${GLNEXUS_TEST_PROJECT_ID}:/
dx build -a --destination :/Attic/travis/GLnexus cli/dxapplet
GIT_REVISION=$(git describe --long --tags --always)
dx set_properties :/Attic/travis/GLnexus "git_revision=${GIT_REVISION}"

# launch on pre-staged test dataset
dx run :/Attic/travis/GLnexus \
    -i gvcf_tar=:/Fake/gvcfs1024.tar -i bed_ranges_to_genotype=:/vcr_b37.bed \
    -i iter_compare=true \
    --name "GLnexus medium tests ${GIT_REVISION}" \
    --folder :/Attic/travis --instance-type mem3_ssd1_x8 --priority normal -y

# We don't wait for completion because typical runs are liable to take longer
# than the travis time limit. You have to look in GLNEXUS_TEST_PROJECT_ID to
# see what happened.
