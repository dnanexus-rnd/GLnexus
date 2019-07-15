# Dockerfile for building GLnexus. The resulting container image runs the unit tests
# by default. It has in its working directory the statically linked glnexus_cli
# executable which can be copied out.
FROM ubuntu:18.04
MAINTAINER DNAnexus
ENV DEBIAN_FRONTEND noninteractive
ARG git_revision=master
ARG build_type=Release

# dependencies
RUN apt-get -qq update && \
     apt-get -qq install -y --no-install-recommends --no-install-suggests \
     curl wget ca-certificates git-core less netbase \
     g++ cmake autoconf make file valgrind \
     libjemalloc-dev libzip-dev libsnappy-dev libbz2-dev zlib1g-dev liblzma-dev libzstd-dev \
     python-pyvcf

# clone GLnexus repo on the desired git revision
WORKDIR /
RUN git clone https://github.com/dnanexus-rnd/GLnexus.git
WORKDIR /GLnexus
RUN git fetch --tags origin && git checkout "$git_revision" && git submodule update --init --recursive

# compile GLnexus
RUN cmake -DCMAKE_BUILD_TYPE=$build_type . && make -j4

# set up default container start to run tests
CMD ctest -V

