#!/bin/bash
set -ex -o pipefail

git_revision=master
if [ -n "$1" ]; then
  git_revision="$1"
fi

# build GLnexus
curl -s https://raw.githubusercontent.com/dnanexus-rnd/GLnexus/${git_revision}/Dockerfile \
    | docker build -t glnexus_tests --build-arg=git_revision=$git_revision -

# run unit tests
docker run --rm glnexus_tests

# copy executables out
rm -rf /tmp/glnexus_docker && mkdir -p /tmp/glnexus_docker
docker run --rm -v /tmp/glnexus_docker:/io glnexus_tests cp glnexus_cli cli/dxapplet/resources/usr/local/bin/bgzip /io

# formulate the tag we'll attach to the final image
tag=$(docker run --rm glnexus_tests bash -c "git fetch --tags origin && git describe --long --always --tags")
tag="quay.io/mlin/glnexus:$tag"

# synthesize a Dockerfile for an image containing glnexus_cli and misc utilities
echo "FROM ubuntu:16.04
MAINTAINER DNAnexus
RUN apt-get -qq update && apt-get -qq install -y libjemalloc-dev numactl bcftools
RUN apt-get clean
COPY glnexus_cli bgzip /usr/local/bin/
" > /tmp/glnexus_docker/Dockerfile

# build image from this synthesized context
docker build --no-cache -t "${tag}-preprecursor" /tmp/glnexus_docker
# flatten the image, to further reduce its deploy size, and set up the runtime ENV/WORKDIR etc.
temp_container_id=$(docker create "${tag}-preprecursor")
docker export "$temp_container_id" | docker import - "${tag}-precursor"
echo "FROM ${tag}-precursor" '
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
WORKDIR /
CMD /bin/bash' | docker build -t "$tag" -

echo "$tag"

