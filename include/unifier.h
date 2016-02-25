#ifndef GLNEXUS_UNIFIER_H
#define GLNEXUS_UNIFIER_H

#include "types.h"

namespace GLnexus {

// unification_config...

/// Compute unified sites from all discovered alleles in some genomic region
Status unified_sites(const unifier_config& cfg,
	                 const discovered_alleles& alleles,
	                 std::vector<unified_site>& ans,
	                 range containing_target = range(-1,-1,-1));

}

#endif
