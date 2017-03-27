#ifndef GLNEXUS_UNIFIER_H
#define GLNEXUS_UNIFIER_H

#include "types.h"

namespace GLnexus {

// unification_config...

/// Compute unified sites from all discovered alleles in some genomic region
/// N = sample count (used to estimate allele frequencies)
// Important I/O notes:
// (1) The input discovered_alleles structure is cleared as a
//     side-effect, to save memory since it can be quite large.
// (2) ans is NOT cleared, only appended to.
Status unified_sites(const unifier_config& cfg,
                     unsigned N,
                     /* const */ discovered_alleles& alleles,
                     std::vector<unified_site>& ans);
}

#endif
