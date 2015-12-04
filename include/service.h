#ifndef GLNEXUS_SERVICE_H
#define GLNEXUS_SERVICE_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include "types.h"
#include "data.h"

namespace GLnexus {

class Service {
    // pImpl idiom
    struct body;
    std::unique_ptr<body> body_;

    Service(BCFData& data);
    Service(const Service&) = delete;

public:
    static Status Start(Metadata& metadata, BCFData& data, std::unique_ptr<Service>& svc);
    ~Service();


    /// Discover all the alleles contained within the given range for a sample
    /// set.

    /// Each allele has a reference range and DNA sequence. Additionally, it's
    /// marked as REF or not, and includes an estimate of the observation
    /// count in the sample set, which may be useful for making an empirical
    /// estimate of the allele frequency. The discovered alleles may be
    /// overlapping; the set is generally used as input to allele
    /// "unification"
    Status discover_alleles(const std::string& sampleset, const range& pos, discovered_alleles& ans);

    /// Discover all the alleles contained within each of the given disjoint
    /// ranges. Uses multithreading such that this may be preferable to
    /// calling discover_alleles repeatedly for many small ranges (e.g. exome
    /// target capture regions). However, attention should be paid to the
    /// anticipated size of the results.
    Status discover_alleles(const std::string& sampleset, const std::vector<range>& ranges,
                            std::vector<discovered_alleles>& ans);

    /// Genotype a set of samples at the given sites, producing a BCF file.
    Status genotype_sites(const genotyper_config& cfg, const std::string& sampleset, const std::vector<unified_site>& sites, const std::string& filename, consolidated_loss& dlosses);
};

}

#endif
