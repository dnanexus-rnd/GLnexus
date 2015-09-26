#!/bin/bash -e

# Build & install valgrind. The Ubuntu 12.04 version from apt is too old.
cd /tmp
curl http://valgrind.org/downloads/valgrind-3.11.0.tar.bz2 | tar xjp
cd valgrind*/
./configure
make -j $(nproc)
sudo make install
