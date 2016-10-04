#ifndef GLNEXUS_COMPARE_QUERY_H
#define GLNEXUS_COMPARE_QUERY_H

#include "BCFKeyValueData.h"

// utilities for comparing two BCF iterators

// generate a random number in the range [0 .. n-1]
int gen_rand_number(int n);

// Generate a number that is one of [0, 1/n, 2/n, 3/n, ... (n-1)/n]
double gen_rand_double(int n);

// Compares the two query iterators, and returns:
//  1: success
//  0: failure, the iterators returned different results
// -1: abort, used too much memory
int compare_query(GLnexus::BCFKeyValueData &data, GLnexus::MetadataCache &cache,
                   const std::string& sampleset, const GLnexus::range& rng);


#endif
