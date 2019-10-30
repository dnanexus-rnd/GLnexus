# GLnexus
**From DNAnexus R&D: scalable gVCF merging and joint variant calling for population sequencing projects.**
(GL, genotype likelihood)

### [bioRxiv preprint](http://dx.doi.org/10.1101/343970)

In our [manuscript](http://dx.doi.org/10.1101/343970) with collaborators at [Regeneron Genetics Center](https://www.regeneron.com/genetics-center) and [Baylor College of Medicine](https://www.hgsc.bcm.edu/), we detail the design of GLnexus and scientific validation using up to 240,000 human exomes and 22,600 genomes. Compared to the DNAnexus cloud-native deployment used for such large projects, this open-source version produces identical scientific results but lacks some of the scalability and production-oriented features.

### [Getting Started](https://github.com/dnanexus-rnd/GLnexus/wiki/Getting-Started)

The [Getting Started](https://github.com/dnanexus-rnd/GLnexus/wiki/Getting-Started) wiki page has a tutorial for first-time users.

### [Prebuilt executables](https://github.com/dnanexus-rnd/GLnexus/releases)

For each tagged revision, the [Releases](https://github.com/dnanexus-rnd/GLnexus/releases) page has a static executable suitable for most Linux x86-64 hosts; just download it and `chmod +x glnexus_cli`. Each release also provides a lightweight Docker image wrapping `glnexus_cli`.

### Build & test

<a href="https://travis-ci.org/dnanexus-rnd/GLnexus"><img src="https://travis-ci.org/dnanexus-rnd/GLnexus.svg?branch=master"/></a> [![Coverage Status](https://coveralls.io/repos/dnanexus-rnd/GLnexus/badge.svg?branch=master&service=github)](https://coveralls.io/github/dnanexus-rnd/GLnexus?branch=master)

The GLnexus build process has a number of dependencies, but produces a standalone, statically-linked executable `glnexus_cli`. The easiest way to build it is to use our Dockerfile to control all the compile-time dependencies, then simply copy the static executable out of the resting Docker container and put it anywhere you like. 

```
# Clone repo
git clone https://github.com/dnanexus-rnd/GLnexus.git
cd GLnexus
git checkout vX.Y.Z  # optional, check out desired revision

# Build GLnexus in docker
docker build --no-cache -t glnexus_tests .

# Run GLnexus unit tests.
docker run --rm glnexus_tests

# Copy the static GLnexus executable to the current working directory.
docker run --rm -v $(pwd):/io glnexus_tests cp glnexus_cli /io

# Run it to see its usage message.
./glnexus_cli
```

**To build GLnexus without Docker**, make sure you have [gcc 5+](http://askubuntu.com/a/581497), [CMake 3.2+](http://askubuntu.com/questions/610291/how-to-install-cmake-3-2-on-ubuntu-14-04), and all the dependencies indicated in the [Dockerfile](https://github.com/dnanexus-rnd/GLnexus/blob/master/Dockerfile). 

Then,

```
git clone https://github.com/dnanexus-rnd/GLnexus.git
cd GLnexus
cmake -Dtest=ON . && make -j$(nproc) && ctest -V
```

You will also find `./glnexus_cli` here.

### Coding conventions

* C++14 - take advantage of [the goodies](http://shop.oreilly.com/product/0636920033707.do)
* Use smart pointers to avoid passing resources needing manual deallocation across function/class boundaries
* Prefer references over pointers when they shouldn't be null nor change ever.
* Avoid exceptions; prefer returning a `Status`, defined early in [types.h](https://github.com/dnanexus-rnd/GLnexus/blob/master/include/types.h)
 * nb the frequently-used convenience macro `S()` defined just below `Status`
* Avoid public constructors with nontrivial bodies; prefer static initializer function returning `Status`
* Avoid elaborate templated class hierarchies

### Libraries used 
* [htslib](https://github.com/samtools/htslib)
* [rocksdb](https://github.com/facebook/rocksdb)
* [yaml-cpp](https://github.com/jbeder/yaml-cpp)
* [Capnproto](https://github.com/sandstorm-io/capnproto)
* [CTPL](https://github.com/vit-vit/CTPL)
* [fcmm](https://github.com/giacomodrago/fcmm)
* [zstd](https://github.com/facebook/zstd)
* [Catch](https://github.com/philsquared/Catch) test framework

### Performance profiling

The [Performance](https://github.com/dnanexus-rnd/GLnexus/wiki/Performance) wiki page has practical advice for deploying GLnexus on a powerful server.

The code has some hooks for performance profiling using
[`perf`](https://en.wikipedia.org/wiki/Perf_(Linux)) and
[FlameGraph](http://www.brendangregg.com/FlameGraphs/cpuflamegraphs.html).

To profile performance within the DNAnexus applet run the applet as
usual plus `-i perf=true`. This produces an output file
```genotype.stacks``` containing sampling observation counts for common call
stacks. To generate an SVG visualization with FlameGraph:

```
git clone https://github.com/brendangregg/FlameGraph
FlameGraph/flamegraph.pl < genotype.stacks > genotype.svg
```
