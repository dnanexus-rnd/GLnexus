#ifndef GLNEXUS_GENOTYPER_H
#define GLNEXUS_GENOTYPER_H

#include "data.h"
#include "types.h"
#include <fstream>
#include <memory>
#include "residuals.h"

namespace GLnexus {

// Genotype a site.
//
// residual_rec: in case there are call losses, generate a YAML formatted record giving
// the context. This is used offline to improve the algorithms.
Status genotype_site(const genotyper_config& cfg, MetadataCache& cache, BCFData& data,
                     const unified_site& site,
                     const std::string& sampleset, const std::vector<std::string>& samples,
                     const bcf_hdr_t* hdr, std::shared_ptr<bcf1_t>& ans,
                     bool residualsFlag,
                     std::shared_ptr<std::string> &residual_rec,
                     std::atomic<bool>* abort = nullptr);

/// A single allele call and metadata; diploid samples each have two calls
struct one_call {
    int32_t allele = bcf_gt_missing; /// or bcf_gt_allele(some_allele)
    NoCallReason RNC = NoCallReason::MissingData;

    one_call() = default;
    one_call(int32_t allele_, NoCallReason RNC_) : allele(allele_), RNC(RNC_) {}

    bool operator==(const one_call& rhs) const noexcept {
        return allele == rhs.allele && RNC == rhs.RNC;
    }

    bool operator < (const one_call& rhs) const noexcept {
        // Ordering of Missing vs Lost calls is arbitrary.
        // TODO: Missing call is set to come after usual calls.
        // Check ordering convention/vcf guidelines
        return (rhs.allele == bcf_gt_missing && allele != bcf_gt_missing) || (allele < rhs.allele) || (allele == rhs.allele && RNC < rhs.RNC);
    }

    bool operator <= (const one_call& rhs) const noexcept {
        return *this < rhs || *this == rhs;
    }

    bool operator > (const one_call& rhs) const noexcept {
        return !(*this <= rhs);
    }
};
} // namespace GLnexus
#endif
