#include "types.h"

using namespace std;

namespace GLnexus {

string bcf_pos_errmsg(bcf1_t *bcf) {
    ostringstream os;
    os << '<' << bcf->rid << '>';
    os << ':' << bcf->pos;
    return os.str();
}

Status range_of_bcf(const bcf_hdr_t* hdr, bcf1_t* bcf, range& ans) {
    if (hdr == nullptr || bcf == nullptr) {
        return Status::Invalid("range_of_bcf null input");
    }

    bcf_info_t *info = bcf_get_info(hdr, bcf, "END");
    if (info) {
        if (info->type != BCF_BT_INT32) {
            return Status::Invalid("range_of_bcf: END field has unexpected type", bcf_pos_errmsg(bcf));
        }
        if (info->len != 1) {
            return Status::Invalid("range_of_bcf: multiple END fields", bcf_pos_errmsg(bcf));
        }
        if (info->v1.i < bcf->pos) {
            return Status::Invalid("range_of_bcf: END < POS", bcf_pos_errmsg(bcf));
        }
        ans.end = info->v1.i;
    } else {
        ans.end = bcf->pos+bcf->rlen;
    }

    ans.rid = bcf->rid;
    ans.beg = bcf->pos;
    return Status::OK();
}

Status range_of_bcf(const bcf_hdr_t* hdr, const std::shared_ptr<bcf1_t>& bcf, range& ans) {
    return range_of_bcf(hdr, bcf.get(), ans);
}

Status range_of_bcf(const std::shared_ptr<const bcf_hdr_t>& hdr,
                    const std::shared_ptr<bcf1_t>& bcf, range& ans) {
    return range_of_bcf(hdr.get(), bcf.get(), ans);
}

// Add src alleles to dest alleles. Identical alleles alleles ale merged,
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

}