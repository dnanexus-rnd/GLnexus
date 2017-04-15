#!/bin/bash -e

# PPAs
sudo add-apt-repository -y ppa:nathan-renniewaldock/ppa # lz4
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test # gcc 4.9
sudo apt-get -qq update

# The lz4 library requires a special repository
sudo apt-get -qq install -y liblz4-dev

# Install GCC 4.9
sudo apt-get -qq install -y gcc-5 g++-5 binutils
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 \
                         --slave /usr/bin/g++ g++ /usr/bin/g++-5 \
                         --slave /usr/bin/gcov gcov /usr/bin/gcov-5

# Install boost for use in yaml-cpp
sudo apt-get -qq install -y libboost1.55-dev
