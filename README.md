# GLnexus
**From DNAnexus R&D: a scalable datastore for population genome sequencing, with on-demand joint genotyping.**
(GL = genotype likelihood)

### Build & run tests

First [install gcc 4.9](http://askubuntu.com/a/581497) and `libjemalloc-dev`. Then:

```
cmake -Dtest=ON . && make && ./unit_tests
```

