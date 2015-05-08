#ifndef GLNEXUS_ALLELES_H
#define GLNEXUS_ALLELES_H

#include "types.h"

namespace GLnexus {

// unification_config...

/// exactly one input allele for any given range must be designated as the reference.
Status unify_alleles(const discovered_alleles& alleles, std::vector<unified_site>& ans);

}

#endif
