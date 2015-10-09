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
    /// Initialize a brand-new database, which SHOULD be empty to begin with.
    /// Contigs are stored and an empty sample set "*" is created.
    static Status InitializeDB(KeyValue::DB* db,
                               const std::vector<std::pair<std::string,size_t> >& contigs,
                               int interval_len = 30000);

    /// Open an existing database
    static Status Open(KeyValue::DB* db, std::unique_ptr<BCFKeyValueData>& ans, size_t cacheBytes = (1 << 30));

    virtual ~BCFKeyValueData();

    // Metadata
    Status contigs(std::vector<std::pair<std::string,size_t> >& ans) const override;
    Status sampleset_samples(const std::string& sampleset,
                             std::shared_ptr<const std::set<std::string> >& ans) const override;
    Status sample_dataset(const std::string& sample, std::string& ans) const override;
    Status all_samples_sampleset(std::string& ans) override;

    // statistics
    std::shared_ptr<StatsRangeQuery> getRangeStats();

    // BCFData
    Status dataset_header(const std::string& dataset,
                              std::shared_ptr<const bcf_hdr_t>& hdr) const override;
    Status dataset_range(const std::string& dataset, const bcf_hdr_t* hdr, const range& pos,
                         std::vector<std::shared_ptr<bcf1_t> >& records) override;

    /// Import a new data set (a gVCF file, possibly containing multiple samples).
    /// The data set name must be unique.
    /// The sample names in the data set (gVCF column names) must be unique.
    /// All samples are immediately added to the sample set "*"
    Status import_gvcf(MetadataCache& metadata, const std::string& dataset,
                       const std::string& filename,
                       std::set<std::string>& samples);
};

}

#endif
