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
/// Sample sets and data sets are immutable.
class Metadata {

public:
    virtual ~Metadata() = default;

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
    virtual Status sample_dataset(const std::string& sample, std::string& ans) const = 0;

    /// Return the name of a sample set representing all samples currently
    /// available. This may either create a new sample set if needed, or
    /// return an existing one if available. As always, the sample set is
    /// immutable: it will not include samples added to the database later
    /// (but one could call all_samples_sampleset again to get a different
    /// sample set including them).
    virtual Status all_samples_sampleset(std::string& ans) = 0;

    /// Return the count of all samples in the database.
    virtual Status sample_count(size_t& ans) const = 0;
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
    Status all_samples_sampleset(std::string& ans) override;
    Status sample_count(size_t& ans) const override;

    const std::vector<std::pair<std::string,size_t> >& contigs() const;
    Status sampleset_datasets(const std::string& sampleset,
                              std::shared_ptr<const std::set<std::string> >& samples,
                              std::shared_ptr<const std::set<std::string>>& datasets) const;
};

/// Iterate over BCF records within some range.
class RangeBCFIterator {
public:
    virtual ~RangeBCFIterator() = default;

    /// Get all the records in one dataset. Returns NotFound at the end of the
    /// iteration.
    virtual Status next(std::string& dataset, std::shared_ptr<const bcf_hdr_t>& hdr,
                        std::vector<std::shared_ptr<bcf1_t>>& records) = 0;
};

/// Abstract interface to stored BCF data sets. The implementation is
/// responsible for any suitable caching.
class BCFData {
public:
    virtual ~BCFData() = default;

    /// Retrieve the BCF header for a data set.
    virtual Status dataset_header(const std::string& dataset,
                                  std::shared_ptr<const bcf_hdr_t>* hdr) = 0;

    /// Retrieve all BCF records in the data set overlapping a range.
    ///
    /// Each record x will already have been "unpacked" with
    /// bcf_unpack(x,BCF_UN_ALL). The records may be shared, so they must not
    /// be mutated. (They aren't declared const because some vcf.h accessor
    /// functions don't take const bcf1_t*)
    ///
    /// predicate: a function that filters BCF records in the range,
    /// only if it returns true, is the record added to [records]. For
    /// example, a typical function return true only if min_alleles is
    /// greater than three.  This is a targeted optimization for
    /// selecting only variant records from gVCF files, as variant
    /// records have at least three alleles (REF, one or more
    /// specified ALTs, and <NON_REF>/<*>) while reference confidence
    /// records have only two alleles. Thus, setting min_alleles=3
    /// should yield variant records only, while min_alleles=0 will
    /// yield all records.
    ///
    /// Note: the BCF record is provided in packed form to the predicate
    /// function. It can unpack it if needed, in which case, it will not be
    /// unpacked again by the dataset_range function.
    ///
    /// The provided header must match the data set, otherwise the behavior is undefined!
    virtual Status dataset_range(const std::string& dataset, const bcf_hdr_t* hdr,
                                 const range& pos, bcf_predicate predicate,
                                 std::vector<std::shared_ptr<bcf1_t>>* records) = 0;

    /// Wrapper for dataset_range which first fetches the appropriate header
    /// (useful if the caller doesn't already have the header in hand)
    virtual Status dataset_range_and_header(const std::string& dataset,
                                            const range& pos, bcf_predicate predicate,
                                            std::shared_ptr<const bcf_hdr_t>* hdr,
                                            std::vector<std::shared_ptr<bcf1_t>>* records);

    /// Get iterators for BCF records overlapping the given range in all
    /// datasets containing at least one sample in the designated sample set.
    //
    /// To facilitate parallelization, the implementation may yield multiple
    /// iterators, each of which will produce a range-based disjoint subset of
    /// the relevant records. Each iterator will yield results for each
    /// relevant data set (possibly yielding zero records in some steps) --
    /// that is, they will all reach their end after the same number of steps.
    /// The iterators together will produce each relevant record exactly once.
    virtual Status sampleset_range(const MetadataCache& metadata, const std::string& sampleset,
                                   const range& pos, bcf_predicate predicate,
                                   std::shared_ptr<const std::set<std::string>>& samples,
                                   std::shared_ptr<const std::set<std::string>>& datasets,
                                   std::vector<std::unique_ptr<RangeBCFIterator>>& iterators);
};

}

#endif
