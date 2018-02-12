# Dockerfile for building GLnexus. The resulting container image runs the unit tests
# by default. It has in its working directory the statically linked glnexus_cli
# executable which can be copied out.
FROM ubuntu:16.04
MAINTAINER DNAnexus
ENV DEBIAN_FRONTEND noninteractive
ARG git_revision=master

# dependencies
RUN apt-get -qq update && \
     apt-get -qq install -y --no-install-recommends --no-install-suggests \
     curl wget ca-certificates git-core less netbase \
     g++ cmake autoconf make file \
     libjemalloc-dev libzip-dev libsnappy-dev libbz2-dev zlib1g-dev liblzma-dev \
     python-pyvcf

# Install Zstandard, which is a rapidly-evolving project -- the versions in the apt
# repositories are too old.
WORKDIR /build
RUN bash -c "tar zxf <(curl -L https://github.com/facebook/zstd/archive/v1.3.2.tar.gz) && make -C zstd-* -j4 && make -C zstd-* install"

# clone GLnexus repo on the desired git revision
RUN git clone https://github.com/dnanexus-rnd/GLnexus.git
WORKDIR /build/GLnexus
RUN git fetch --tags origin && git checkout "$git_revision" && git submodule update --init --recursive

# compile GLnexus
RUN cmake -Dtest=ON -DCMAKE_BUILD_TYPE=Release . && make -j4

# set up default container start to run tests
CMD ./unit_tests

