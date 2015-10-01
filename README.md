# GLnexus
**From DNAnexus R&D: a scalable datastore for population genome sequencing, with on-demand joint genotyping.**
(GL, genotype likelihood)

This is an early-stage R&D project we're developing openly. The code doesn't yet do anything useful! There's a [wiki project roadmap](https://github.com/dnanexus-rnd/GLnexus/wiki), which should be read in the spirit of "plans are worthless, but planning is insdispensable."

### Build & run tests

<a href="https://travis-ci.org/dnanexus-rnd/GLnexus"><img src="https://travis-ci.org/dnanexus-rnd/GLnexus.svg?branch=master"/></a> [![Coverage Status](https://coveralls.io/repos/dnanexus-rnd/GLnexus/badge.svg?branch=master&service=github)](https://coveralls.io/github/dnanexus-rnd/GLnexus?branch=master)

First [install gcc 4.9](http://askubuntu.com/a/581497), `cmake`, `libjemalloc-dev`, `libboost-dev`, `libzip-dev`, `libsnappy-dev`, `liblz4-dev`, `libbz2-dev`. Then:

```
cmake -Dtest=ON . && make && ./unit_tests
```

Other dependencies (should be set up automatically by CMake):
* [htslib](https://github.com/samtools/htslib)
* [rocksdb](https://github.com/facebook/rocksdb)
* [yaml-cpp](https://github.com/jbeder/yaml-cpp)
* [Catch](https://github.com/philsquared/Catch) test framework

### Developer documentation

Evolving developer documentation can be found on the [project github page](http://dnanexus-rnd.github.io/GLnexus/index.html).

### Coding conventions

* C++14 - take advantage of [the goodies](http://shop.oreilly.com/product/0636920033707.do)
* Use smart pointers to avoid passing resources needing manual deallocation across function/class boundaries
* Prefer references over pointers when they shouldn't be null nor change ever.
* Avoid exceptions; prefer returning a `Status`, defined early in [types.h](https://github.com/dnanexus-rnd/GLnexus/blob/master/include/types.h)
 * nb the frequently-used convenience macro `S()` defined just below `Status`
* Avoid public constructors with nontrivial bodies; prefer static initializer function returning `Status`
* Avoid elaborate templated class hierarchies


### How to do performance profiling

In order to do performance profiling on the code, you need
to compile it with special c++ flags. This can be done as follows:
```
cmake -Dtest=ON -DCMAKE_BUILD_TYPE=Profile .
```

Running the glnexus_cli applet produces several outputs, one of which
is a file called ```genotype.stacks```. It contains a list of common
call stacks executed by the glnexus code, download it to your
local machine. To manipulate this file, please install the FlameGraphs
package. Here is an example for generating an SVG graph from the run results:

```
git clone https://github.com/brendangregg/FlameGraph
grep glnexus genotype.stacks > X
FlameGraph/flamegraph.pl --minwidth 0.5 X > genotype.svg
```
