#include "types.h"

using namespace std;

namespace GLnexus {

// Add src alleles to dest alleles. Identical alleles alleles ale merged,
// using the sum of their observation_counts
Status merge_discovered_alleles(const discovered_alleles& src, discovered_alleles& dest) {
    for (auto& dsal : src) {
    	UNPAIR(dsal,allele,ai)
    	UNPAIR(ai,is_ref,obs_count)
    	auto p = dest.find(allele);
        if (p == dest.end()) {
            dest[allele] = ai;
        } else {
        	if (is_ref != p->second.first) {
        		return Status::Invalid("allele appears as both REF and ALT", allele.dna + "@" + allele.pos.str());
        	}
            dest[allele] = make_pair(is_ref,obs_count+p->second.second);
        }
    }

    return Status::OK();
}

}