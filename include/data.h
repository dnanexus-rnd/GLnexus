#ifndef GLNEXUS_DATA_H
#define GLNEXUS_DATA_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <set>
#include <vcf.h>
#include "types.h"

namespace GLnexus {

/// Abstract interface to "metadata" -- namely contigs and sample sets.
///
/// Definitions:
/// Each sample (individual whose GLs are stored) has a unique string identifier.
/// A sample set is simply a set of sample identifiers.
/// A data set contains all the allele/genotype/likelihood data for one or more
/// samples -- just like a [V/B]CF file.
/// All the GL data for each sample resides in exactly one data set.
/// Sample sets and data sets are immutable once created/sealed, except
/// possibly for a special "*" sample set representing all available samples.
class Metadata {

public:
    /// Get the reference contigs.
    ///
    /// The indices of the vector are the "rid" used in range()
    virtual Status contigs(std::vector<std::pair<std::string,size_t> >& ans) const = 0;

    /// List the samples in a sample set.
    ///
    /// The resulting data structure may be shared, so the strings must not be
    /// mutated. They aren't declared const because...C++
    /// http://stackoverflow.com/a/21365478
    virtual Status sampleset_samples(const std::string& sampleset,
                                     std::shared_ptr<const std::set<std::string> >& ans) const = 0;

    /// Find the data set containing the sample.
    ///
    /// The data set may contain other samples.
    /// 
    /// The resulting data structure may be shared, so the strings must not be
    /// mutated. They aren't declared const because...C++
    /// http://stackoverflow.com/a/21365478
    virtual Status sample_dataset(const std::string& sample, std::string& ans) const = 0;
};


/// Wraps any Metadata implementation to provide in-memory caching/indexing of
/// the immutable relationships
class MetadataCache : public Metadata {
    struct body;
    std::unique_ptr<body> body_;

    MetadataCache();
    MetadataCache(const MetadataCache&) = delete;

public:
    static Status Start(Metadata& inner, std::unique_ptr<MetadataCache>& ptr);
    virtual ~MetadataCache();

    Status contigs(std::vector<std::pair<std::string,size_t> >& ans) const override;
    Status sampleset_samples(const std::string& sampleset,
                             std::shared_ptr<const std::set<std::string> >& ans) const override;
    Status sample_dataset(const std::string& sample, std::string& ans) const override;

    const std::vector<std::pair<std::string,size_t> >& contigs() const;
    Status sampleset_datasets(const std::string& sampleset,
                              std::shared_ptr<const std::set<std::string> >& samples,
                              std::shared_ptr<const std::set<std::string>>& datasets) const;
};

/// Abstract interface to stored BCF data sets. The implementation is
/// responsible for any suitable caching.
class BCFData {
public: 
    /// Retrieve the BCF header for a data set.
    virtual Status dataset_header(const std::string& dataset,
                                  std::shared_ptr<const bcf_hdr_t>& hdr) const = 0;

    /// Retrieve all BCF records in the data set overlapping a range.
    ///
    /// Each record x will already have been "unpacked" with
    /// bcf_unpack(x,BCF_UN_ALL). The records may be shared, so they must not
    /// be mutated. (They aren't declared const because some vcf.h accessor
    /// functions don't take const bcf1_t*)
    ///
    /// The provided header must match the data set, otherwise the behavior is undefined!
    virtual Status dataset_range(const std::string& dataset, const bcf_hdr_t* hdr, const range& pos,
                                 std::vector<std::shared_ptr<bcf1_t> >& records) const = 0;

    /// Wrapper for dataset_range which first fetches the appropriate header
    /// (useful if the caller doesn't already have the header in hand)
    virtual Status dataset_range_and_header(const std::string& dataset, const range& pos,
                                            std::shared_ptr<const bcf_hdr_t>& hdr,
                                            std::vector<std::shared_ptr<bcf1_t> >& records) const;

};

}

#endif
