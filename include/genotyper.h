#ifndef GLNEXUS_GENOTYPER_H
#define GLNEXUS_GENOTYPER_H

#include "data.h"
#include "types.h"
#include <memory>
#include "genotyper_config.h"

namespace GLnexus {

Status genotype_site(const genotyper_config& cfg, BCFData& data, const unified_site& site,
                     const std::vector<std::string>& samples, const std::set<std::string>& datasets,
                     const bcf_hdr_t* hdr, std::shared_ptr<bcf1_t>& ans, consolidated_loss& losses_for_site);

// LossTracker handles the low-level housekeeping of loss accounting for a
// single unified_site and a single sample. Computation in this class assumes a
// diploid genome.
// The get() function returns the computed loss as a loss_stats for higher
// level record keeping.
class LossTracker {

public:
    // Simple wrapper struct to store information about an original call
    // for a unified site, used bfor computation of loss
    struct orig_call {
        orig_call(range pos_, bool is_gvcf_) : pos(pos_), is_gvcf(is_gvcf_) {}

        range pos;
        bool is_gvcf;

        bool operator==(const orig_call& rhs) const noexcept { return pos == rhs.pos && is_gvcf == rhs.is_gvcf; }
        bool operator<(const orig_call& rhs) const noexcept { return pos < rhs.pos; }
        bool operator<=(const orig_call& rhs) const noexcept { return pos <= rhs.pos; }
    };

    // Constructor
    LossTracker(const range rng_) noexcept : rng(rng_) {}

    Status add_call_for_site(const range call, int n_calls, bool is_gvcf) noexcept;
    Status finalize_loss_for_site(int n_no_calls) noexcept;
    Status get(loss_stats& ans) const noexcept;

private:
    // Range of joint-called unified_site being considered
    range rng;

    // Original calls (identified by effective range within site) and
    // count of calls. Calls which are different in the bcf record but
    // share the same effective range within the site will be collapsed
    // into the same key.
    std::map<orig_call, int> orig_calls_for_site;

    int n_calls_total=0, n_bp_total=0;
    int n_gvcf_calls_total=0, n_gvcf_bp_total=0;
    int n_calls_lost=0, n_bp_lost=0;
    int n_gvcf_calls_lost=0, n_gvcf_bp_lost=0;
    int n_no_calls_total = 0;

    bool is_finalized = false;
};

using LossTrackers = std::vector<LossTracker>;

} // namespace GLnexus
#endif
