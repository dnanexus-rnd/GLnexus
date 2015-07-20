#ifndef GLNEXUS_BCFKEYVALUEDATA_H
#define GLNEXUS_BCFKEYVALUEDATA_H
#include "data.h"
#include "KeyValue.h"

namespace GLnexus {

/// An implementation of the Data interface that stores sample sets and BCF
/// records in a given key-value database. One imported gVCF file (potentially
/// with multiple samples) becomes a data set. The key schema permits
/// efficient retrieval by genomic range across the datasets.
class BCFKeyValueData : public Data {
    // pImpl idiom
    struct body;
    std::unique_ptr<body> body_;

    BCFKeyValueData();
public:
    /// Initialize a brand-new database, which SHOULD be empty.
    static Status InitializeDB(KeyValue::DB* db, const std::vector<std::pair<std::string,size_t> >& contigs);

    /// Open an existing database
    static Status Open(KeyValue::DB* db, std::unique_ptr<BCFKeyValueData>& ans);

    virtual ~BCFKeyValueData();

    // Data interface
    Status contigs(std::vector<std::pair<std::string,size_t> >& ans) const override;
    Status sampleset_samples(const std::string& sampleset,
                             std::shared_ptr<const std::set<std::string> >& ans) const override;
    Status sample_dataset(const std::string& sample, std::string& ans) const override;
    
    Status dataset_bcf_header(const std::string& dataset,
                              std::shared_ptr<const bcf_hdr_t>& hdr) const override;
    Status dataset_bcf(const std::string& dataset, const range& pos,
                       std::shared_ptr<const bcf_hdr_t>& hdr,
                       std::vector<std::shared_ptr<bcf1_t> >& records) const override;

    /// Import a new dataset
    Status import_gvcf(const DataCache* cache, const std::string& dataset, const std::string& filename);
};

}

#endif
