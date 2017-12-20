#ifndef GLNEXUS_UNIFIER_H
#define GLNEXUS_UNIFIER_H

#include "types.h"

namespace GLnexus {

// unification_config...
struct unifier_stats {
    // # ALT alleles represented in idiomatic sites.
    size_t unified_alleles = 0;

    // # ALT alleles passing quality thresholds but 'lost' due to not unifying
    // cleanly with other overlapping alleles. If
    // unifier_config::monoallelic_sites_for_lost_alleles is true, then these
    // are included as monoallelic sites (but still NOT counted in
    // unified_alleles.)
    size_t lost_alleles = 0;

    // # alleles discarded due to failing quality thresholds
    size_t filtered_alleles = 0;

    unifier_stats& operator+=(const unifier_stats& rhs) {
        unified_alleles += rhs.unified_alleles;
        lost_alleles += rhs.lost_alleles;
        filtered_alleles += rhs.filtered_alleles;
        return *this;
    }
};

/// Compute unified sites from all discovered alleles in some genomic region
/// N = sample count (used to estimate allele frequencies)
/// target_ranges: If nonempty, only emit sites overlapping these ranges.
// Important I/O notes:
// (1) The input discovered_alleles structure is cleared as a
//     side-effect, to save memory since it can be quite large.
// (2) ans is NOT cleared, only appended to.
Status unified_sites(const unifier_config& cfg,
                     unsigned N,
                     /* const */ discovered_alleles& alleles,
                     const std::set<range>& target_ranges,
                     std::vector<unified_site>& ans,
                     unifier_stats& stats);

// Find which range overlaps [pos]. The ranges are assumed to be non-overlapping.
// (exposed for unit testing)
Status find_target_range(const std::set<range> &ranges, const range &pos, range &ans);

}

#endif
