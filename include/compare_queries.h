#ifndef GLNEXUS_COMPARE_QUERIES_H
#define GLNEXUS_COMPARE_QUERIES_H

#include "BCFKeyValueData.h"

namespace GLnexus {
namespace compare_queries {
// utilities for comparing two BCF iterators.
//
// We want to make sure that our efficient, multi-threaded, scanning code
// is correct. We do this by comparing it to a much simpler iterator.

// generate a random number in the range [0 .. n-1]
int gen_rand_number(int n);

// Generate a number that is one of [0, 1/n, 2/n, 3/n, ... (n-1)/n]
double gen_rand_double(int n);

// Compares the two query iterators, and returns:
//  1: success
//  0: failure, the iterators returned different results
// -1: abort, used too much memory
int compare_query(BCFKeyValueData &data, MetadataCache &cache,
                  const std::string& sampleset, const range& rng);

Status compare_n_queries(int n_iter,
                         BCFKeyValueData &data,
                         MetadataCache &metadata,
                         const std::string& sampleset);
}}

#endif
