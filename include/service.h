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

struct service_config {
    size_t threads = 0;
};

class Service {
    // pImpl idiom
    struct body;
    std::unique_ptr<body> body_;

    Service(const service_config& cfg, BCFData& data);
    Service(const Service&) = delete;

public:
    static Status Start(const service_config& cfg, Metadata& metadata, BCFData& data,
                        std::unique_ptr<Service>& svc);
    ~Service();


    /// Discover all the alleles contained within the given range for a sample
    /// set.

    /// Each allele has a reference range and DNA sequence. Additionally, it's
    /// marked as REF or not, and includes an estimate of the observation
    /// count in the sample set, which may be useful for making an empirical
    /// estimate of the allele frequency. The discovered alleles may be
    /// overlapping; the set is generally used as input to allele
    /// "unification"
    ///
    /// Sets N = sample count (size of sampleset)
    ///
    /// If the abort flag is non-null, then from time to time the algorithm
    /// checks its contents; if it finds the flag set (presumably by some
    /// other thread), then it will cleanly shut down and return status
    /// ABORTED.
    Status discover_alleles(const std::string& sampleset, const range& pos,
                            unsigned& N, discovered_alleles& ans,
                            std::atomic<bool>* abort = nullptr);

    /// Discover all the alleles contained within each of the given disjoint
    /// ranges. Uses multithreading such that this may be preferable to
    /// calling discover_alleles repeatedly for many small ranges (e.g. exome
    /// target capture regions). However, attention should be paid to the
    /// anticipated size of the results.
    Status discover_alleles(const std::string& sampleset, const std::vector<range>& ranges,
                            unsigned& N, std::vector<discovered_alleles>& ans,
                            std::atomic<bool>* abort = nullptr);

    /// Genotype a set of samples at the given sites, producing a BCF file.
    Status genotype_sites(const genotyper_config& cfg, const std::string& sampleset,
                          const std::vector<unified_site>& sites,
                          const std::string& filename,
                          std::atomic<bool>* abort = nullptr);

    // Report cumulative time (milliseconds) worker threads in the above
    // operations have spent 'stalled' waiting on single-threaded processing
    // steps (e.g. output serialization)
    unsigned threads_stalled_ms() const;
};

}

#endif
