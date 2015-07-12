# GLnexus
**From DNAnexus R&D: a scalable datastore for population genome sequencing, with on-demand joint genotyping.**
(GL, genotype likelihood)

This is an early-stage R&D project we're developing openly. The code doesn't yet do anything useful! There's a [wiki project roadmap](https://github.com/dnanexus-rnd/GLnexus/wiki), which should be read in the spirit of "plans are worthless, but planning is insdispensable."

### Build & run tests

First [install gcc 4.9](http://askubuntu.com/a/581497), `libjemalloc-dev`, `libboost-dev`. Then:

```
cmake -Dtest=ON . && make && ./unit_tests
```

### Coding conventions

* C++14, CMake
* Dependencies (should be set up automatically by CMake)
 * [htslib](https://github.com/samtools/htslib)
 * [rocksdb](https://github.com/facebook/rocksdb)
 * [yaml-cpp](https://github.com/jbeder/yaml-cpp)
 * [Catch](https://github.com/philsquared/Catch) test framework
* Avoid exceptions; prefer returning a `Status`, defined early in [types.h](https://github.com/dnanexus-rnd/GLnexus/blob/master/include/types.h)
 * nb the frequently-used convenience macro `S()` defined just below `Status`
* Avoid elaborate templated class hierarchies
* Prefer static initializer function returning `Status` over public constructor for any nontrivial body
