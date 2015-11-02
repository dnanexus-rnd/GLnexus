#!/bin/bash -e

# PPAs
sudo add-apt-repository -y ppa:libreoffice/ppa # doxygen
sudo add-apt-repository -y ppa:nathan-renniewaldock/ppa # lz4
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test # gcc 4.9
sudo apt-get -qq update

# Install current version of doxygen
sudo apt-get -qq install -y doxygen

# The lz4 library requires a special repository
sudo apt-get -qq install -y liblz4-dev

# install updated version of cmake
wget http://www.cmake.org/files/v3.3/cmake-3.3.0-rc3-Linux-x86_64.sh
sh cmake-3.3.0-rc3-Linux-x86_64.sh --prefix=$HOME --exclude-subdir

# Install GCC 4.9
sudo apt-get -qq install -y gcc-4.9 g++-4.9 binutils
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 60 \
                         --slave /usr/bin/g++ g++ /usr/bin/g++-4.9 \
                         --slave /usr/bin/gcov gcov /usr/bin/gcov-4.9

# Install boost for use in yaml-cpp
sudo apt-get -qq install -y libboost1.55-dev

# Install cpp-coveralls
sudo pip install cpp-coveralls
