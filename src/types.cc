#include "types.h"

using namespace std;

namespace GLnexus {

// Add src alleles to dest alleles. Identical alleles alleles are merged,
// using the sum of their observation_counts
Status merge_discovered_alleles(const discovered_alleles& src, discovered_alleles& dest) {
    for (auto& dsal : src) {
        UNPAIR(dsal,allele,ai)
        auto p = dest.find(allele);
        if (p == dest.end()) {
            dest[allele] = ai;
        } else {
            if (ai.is_ref != p->second.is_ref) {
                return Status::Invalid("allele appears as both REF and ALT", allele.dna + "@" + allele.pos.str());
            }
            p->second.observation_count += ai.observation_count;
        }
    }

    return Status::OK();
}


Status merge_loss_stats(const consolidated_loss& src, consolidated_loss& dest) {

    Status s;
    for (auto& sample_stats : src) {
        UNPAIR(sample_stats, sample, src_stats)
        auto target = dest.find(sample);
        if (target == dest.end()) {
            // sample not found in destination
            dest.insert(make_pair(sample, src_stats));
        } else {
            S(target->second.merge_loss(src_stats));
        }
    }
    return Status::OK();
}
}
