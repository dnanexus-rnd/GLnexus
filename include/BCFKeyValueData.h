#ifndef GLNEXUS_BCFKEYVALUEDATA_H
#define GLNEXUS_BCFKEYVALUEDATA_H
#include "data.h"
#include "KeyValue.h"
#include "BCFSerialize.h"

namespace GLnexus {

struct BCFKeyValueData_body;

/// Implements the Metadata and BCFData interfaces with everything stored in a
/// given key-value database. One imported gVCF file (potentially with
/// multiple samples) becomes a data set. The key schema permits efficient
/// retrieval by genomic range across the datasets.
class BCFKeyValueData : public Metadata, public BCFData {
private:
    // pImpl idiom
    std::unique_ptr<BCFKeyValueData_body> body_;

    BCFKeyValueData();
    BCFKeyValueData(const BCFKeyValueData&) = delete;

public:
    static const int default_bucket_size = 3000;

    /// Initialize a brand-new database, which SHOULD be empty to begin with.
    /// Contigs are stored and an empty sample set "*" is created.
    static Status InitializeDB(KeyValue::DB* db,
                               const std::vector<std::pair<std::string,size_t> >& contigs,
                               int interval_len = default_bucket_size);

    /// Open an existing database
    static Status Open(KeyValue::DB* db, std::unique_ptr<BCFKeyValueData>& ans);

    virtual ~BCFKeyValueData();

    // Metadata
    Status contigs(std::vector<std::pair<std::string,size_t> >& ans) const override;
    Status sampleset_samples(const std::string& sampleset,
                             std::shared_ptr<const std::set<std::string> >& ans) const override;
    Status sample_dataset(const std::string& sample, std::string& ans) const override;
    Status all_samples_sampleset(std::string& ans) override;
    Status sample_count(size_t& ans) const override;

    Status new_sampleset(MetadataCache& metadata, const std::string& sampleset,
                         const std::set<std::string>& samples);

    // statistics
    std::shared_ptr<StatsRangeQuery> getRangeStats();

    // BCFData
    Status dataset_header(const std::string& dataset,
                              std::shared_ptr<const bcf_hdr_t>& hdr) const override;
    Status dataset_range(const std::string& dataset, const bcf_hdr_t* hdr,
                         const range& pos, bcf_predicate predicate,
                         std::vector<std::shared_ptr<bcf1_t> >& records) override;

    Status sampleset_range(const MetadataCache& metadata, const std::string& sampleset,
                           const range& pos, bcf_predicate predicate,
                           std::shared_ptr<const std::set<std::string>>& samples,
                           std::shared_ptr<const std::set<std::string>>& datasets,
                           std::vector<std::unique_ptr<RangeBCFIterator>>& iterators) override;

    // Provide a way to call the non-optimized base implementation of
    // sampleset_range. Mostly for unit testing.
    Status sampleset_range_base(const MetadataCache& metadata, const std::string& sampleset,
                                const range& pos, bcf_predicate predicate,
                                std::shared_ptr<const std::set<std::string>>& samples,
                                std::shared_ptr<const std::set<std::string>>& datasets,
                                std::vector<std::unique_ptr<RangeBCFIterator>>& iterators);

    struct import_result {
        std::set<std::string> samples;
        uint64_t records = 0;     // total # BCF records
        uint64_t max_records = 0; // max # records in any bucket
        size_t bytes = 0;             // total BCF bytes
        size_t max_bytes = 0;         // max bytes in any bucket
        uint64_t buckets = 0;     // # buckets
        uint64_t duplicate_records = 0; // # of records duplicated in multiple buckets
        uint64_t skipped_records = 0; // # of records skipped in source gVCF for various caller-specific reasons

        import_result& add_bucket(uint64_t bucket_records, size_t bucket_bytes, uint64_t duplicates) {
            records += bucket_records;
            max_records = std::max(max_records, bucket_records);
            bytes += bucket_bytes;
            max_bytes = std::max(max_bytes, bucket_bytes);
            buckets++;
            duplicate_records += duplicates;
            return *this;
        }

        import_result& operator+=(const import_result& rhs) {
            #ifndef NDEBUG
            size_t sz0 = samples.size();
            #endif
            samples.insert(rhs.samples.begin(), rhs.samples.end());
            assert(samples.size() == sz0+rhs.samples.size());
            records += rhs.records;
            max_records = std::max(max_records, rhs.max_records);
            bytes += rhs.bytes;
            max_bytes = std::max(max_bytes, rhs.max_bytes);
            buckets += rhs.buckets;
            duplicate_records += rhs.duplicate_records;
            skipped_records += rhs.skipped_records;
            return *this;
        }
    };

    /// Import a new data set (a gVCF file, possibly containing multiple samples).
    /// The data set name must be unique.
    /// The sample names in the data set (gVCF column names) must be unique.
    /// All samples are immediately added to the sample set "*"
    /// If range_filter is nonempty, then import only records overlapping one
    /// of those ranges.
    Status import_gvcf(MetadataCache& metadata, const std::string& dataset,
                       const std::string& filename,
                       const std::set<range>& range_filter,
                       import_result& rslt);

    Status import_gvcf(MetadataCache& metadata, const std::string& dataset,
                       const std::string& filename,
                       std::set<std::string>& samples_imported) {
        import_result rslt;
        Status s = import_gvcf(metadata, dataset, filename, {}, rslt);
        samples_imported = move(rslt.samples);
        return s;
    }
};

/// Get the bucket key prefix length for the bcf collection. This is used with
/// e.g. the RocksDB prefix mode, or DynamoDB hash-range key.
size_t BCFKeyValueDataPrefixLength();

}

#endif
