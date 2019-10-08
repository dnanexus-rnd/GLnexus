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
//
// May set ans to nullptr if the site ends up with all ALT alleles trimmed.
Status genotype_site(const genotyper_config& cfg, MetadataCache& cache, BCFData& data,
                     const unified_site& site,
                     const std::string& sampleset, const std::vector<std::string>& samples,
                     const bcf_hdr_t* hdr, std::shared_ptr<bcf1_t>& ans,
                     bool residualsFlag,
                     std::shared_ptr<std::string> &residual_rec,
                     std::atomic<bool>* abort = nullptr);

// Reasons for emitting a non-call (.), encoded in the RNC FORMAT field in the
// output VCF
enum class NoCallReason {
    N_A,                  /// not applicable (the genotype *is* called)
    MissingData,          /// no gVCF coverage at all
    PartialData,          /// partial gVCF coverage
    InsufficientDepth,    /// insufficient depth of coverage
    LostDeletion,         /// unrepresentable overlapping deletion
    LostAllele,           /// unrepresentable allele (other than overlapping deletion)
    UnphasedVariants,     /// site spans multiple unphased variants
    OverlappingVariants,  /// site spans multiple variants which overlap each other
    MonoallelicSite,      /// site is monoallelic; no assertion about the presence of either ref or alt allele
    InputNonCalled,       /// the relevant input gVCF record is itself non-called
};

/// A single allele call and metadata; diploid samples each have two calls
struct one_call {
    int32_t allele = bcf_gt_missing; /// or bcf_gt_allele(some_allele)
    NoCallReason RNC = NoCallReason::MissingData;
    bool half_call = false;

    one_call() = default;
    one_call(int32_t allele_, NoCallReason RNC_) : allele(allele_), RNC(RNC_) {}

    bool operator==(const one_call& rhs) const noexcept {
        return allele == rhs.allele && RNC == rhs.RNC && half_call == rhs.half_call;
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

// exposed for unit testing
struct bcf1_t_plus {
    std::shared_ptr<bcf1_t> p;

    // is_gvcf_ref_record(p)
    bool is_ref = true;

    // GT vector
    htsvecbox<int> gt;

    // mapping from bcf1_t's alleles onto unified alleles
    std::vector<int> allele_mapping;

    // indicating whether each allele is a deletion wrt the reference
    std::vector<bool> deletion_allele;

    // True if the original record was 'haploid' (e.g. GT=0 or GT=1) as seen
    // with some gVCF callers like Strelka2.
    // In preprocessing, we rewrite gt to loook like ./0 or ./1 for simpler
    // treatment subsequently, but note the original form here.
    bool was_haploid = false;
};
Status preprocess_record(const unified_site& site, const bcf_hdr_t* hdr, const std::shared_ptr<bcf1_t>& record, bcf1_t_plus& ans);
Status revise_genotypes(const genotyper_config& cfg, const unified_site& us, const std::map<int, int>& sample_mapping,
                        const bcf_hdr_t* hdr, bcf1_t_plus& vr);

} // namespace GLnexus
#endif
