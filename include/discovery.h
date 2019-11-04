#pragma once
#include "data.h"
#include "types.h"

// Algorithms used in allele discovery.

namespace GLnexus {

// Discover alleles from a RangeBCFIterator. Records not contained within pos will be ignored.
Status discover_alleles_from_iterator(const std::set<std::string>& samples,
                                      const range& pos,
                                      RangeBCFIterator& iterator,
                                      discovered_alleles& dsals);

// verify that the discovered_alleles has a REF allele for each ALT allele, and that the REF
// alleles are all consistent with each other.
Status discovered_alleles_refcheck(const discovered_alleles& als,
                                   const std::vector<std::pair<std::string,size_t>>& contigs);

} // namespace GLnexus
