#!/bin/bash

# Install current version of doxygen
sudo add-apt-repository -y ppa:libreoffice/ppa
sudo apt-get update
sudo apt-get install doxygen

# The lz4 library requires a special repository
sudo add-apt-repository -y ppa:nathan-renniewaldock/ppa
sudo apt-get update
sudo apt-get install liblz4-dev

# install updated version of cmake
wget http://www.cmake.org/files/v3.3/cmake-3.3.0-rc3-Linux-x86_64.sh
sh cmake-3.3.0-rc3-Linux-x86_64.sh --prefix=$HOME --exclude-subdir

# Install GCC 4.9
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gcc-4.9 g++-4.9
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.9

# Install somewhat current version of boost for use in yaml-cpp
sudo add-apt-repository -y ppa:boost-latest/ppa
sudo apt-get update
sudo apt-get install boost1.55

