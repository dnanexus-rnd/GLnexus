# GLnexus
**From DNAnexus R&D: a scalable datastore for gVCF merging and joint variant calling in population genome sequencing projects.**
(GL, genotype likelihood)

This is an early-stage R&D project we're developing openly. The code doesn't yet do anything useful! There's a [wiki project roadmap](https://github.com/dnanexus-rnd/GLnexus/wiki), which should be read in the spirit of "plans are worthless, but planning is indispensable."

<a href="https://travis-ci.org/dnanexus-rnd/GLnexus"><img src="https://travis-ci.org/dnanexus-rnd/GLnexus.svg?branch=master"/></a> [![Coverage Status](https://coveralls.io/repos/dnanexus-rnd/GLnexus/badge.svg?branch=master&service=github)](https://coveralls.io/github/dnanexus-rnd/GLnexus?branch=master)

### Build & test (Linux)

The GLnexus build process has a number of dependencies, but produces a standalone, statically-linked executable `glnexus_cli`. The easiest way to build it is to use our Dockerfile to control all the compile-time dependencies, then simply copy the static executable out of the resting Docker container and put it anywhere you like. 

```
# Build GLnexus using its Dockerfile.
# You can set a specific git revision by adding --build-arg=git_revision=xxxx
curl -s https://raw.githubusercontent.com/dnanexus-rnd/GLnexus/master/Dockerfile | docker build -t glnexus_tests -

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
git clone --recursive https://github.com/dnanexus-rnd/GLnexus.git
cd GLnexus
cmake -Dtest=ON . && make -j4 && ./unit_tests
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
