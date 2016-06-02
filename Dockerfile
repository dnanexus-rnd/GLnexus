FROM ubuntu:14.04
MAINTAINER DNAnexus R&D

ENV DEBIAN_FRONTEND noninteractive

# Install all dependencies required for building GLnexus
RUN apt-get update && apt-get install -y software-properties-common && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && apt-get update && \ 
    apt-get install -y gcc-4.9 g++-4.9 && \
    apt-get install -y cmake libjemalloc-dev libboost-dev libzip-dev libsnappy-dev liblz4-dev libbz2-dev && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.9 && \
    echo gcc --version && \
    echo g++ --version && \
    apt-get install -y git

# Build GLnexus
RUN git clone https://github.com/dnanexus-rnd/GLnexus.git && \
    cd GLnexus && cmake -Dtest=ON . && make -j

# Set up environment for unit tests and run them
RUN apt-get install -y python python-pip && \
    pip install PyVCF && \
    cd GLnexus && ./unit_tests

# Set up environment for convenient use in container
