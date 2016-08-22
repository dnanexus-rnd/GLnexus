#ifndef GLNEXUS_UNIFIER_H
#define GLNEXUS_UNIFIER_H

#include "types.h"

namespace GLnexus {

// unification_config...

/// Compute unified sites from all discovered alleles in some genomic region
/// N = sample count (used to estimate allele frequencies)
Status unified_sites(const unifier_config& cfg,
                     unsigned N, const discovered_alleles& alleles,
                     std::vector<unified_site>& ans);
}

#endif
