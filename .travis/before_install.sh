#!/bin/bash -e

# Install gcc 5
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get -qq update
sudo apt-get -qq install -y gcc-6 g++-6 binutils
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 60 \
                         --slave /usr/bin/g++ g++ /usr/bin/g++-6 \
                         --slave /usr/bin/gcov gcov /usr/bin/gcov-6

# Install boost for use in yaml-cpp
sudo apt-get -qq install -y libboost1.55-dev

# Install zstd
tar zxf <(curl -L https://github.com/facebook/zstd/archive/v1.3.4.tar.gz)
make -C zstd-* -j8
sudo make -C zstd-* install
