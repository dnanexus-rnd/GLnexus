# GLnexus
**From DNAnexus R&D: a scalable datastore for population genome sequencing, with on-demand joint genotyping.**
(GL = genotype likelihood)

### Build & run tests

First [install gcc 4.9](http://askubuntu.com/a/581497), `libjemalloc-dev`, `libboost-dev`. Then:

```
cmake -Dtest=ON . && make && ./unit_tests
```

### Coding conventions

* C++15, CMake
* Dependencies (should be set up automatically by CMake)
 * [Catch](https://github.com/philsquared/Catch) test framework
 * [htslib](https://github.com/samtools/htslib)
 * [rocksdb](https://github.com/facebook/rocksdb)
 * [yaml-cpp](https://github.com/jbeder/yaml-cpp)
* Avoid using exceptions
* Instead of exceptions, prefer returning a `Status` object, which is defined early in [types.h](https://github.com/dnanexus-rnd/GLnexus/blob/master/include/types.h)
 * nb the frequently-used convenience macro `S()` defined just below `Status`
* Prefer static initializer function returning `Status` over public constructor for any nontrivial body
