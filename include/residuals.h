#ifndef GLNEXUS_RESIDUALS_H
#define GLNEXUS_RESIDUALS_H

#include <fstream>
#include <memory>
#include "data.h"
#include "types.h"

namespace GLnexus {

// Dataset and the records it has for a particular site
struct DatasetSiteInfo {
    std::string name;
    std::shared_ptr<const bcf_hdr_t> header;
    std::vector<std::shared_ptr<bcf1_t>> records;
};

// Create a YAML node describing a loss. The node is formatted as
// a string for simplicity. The [sites] variable is a list of
// sites, and records in them, to print out.
Status residuals_gen_record(const unified_site& site,
                            const bcf_hdr_t *gl_hdr,
                            const bcf1_t *gl_call,
                            const std::vector<DatasetSiteInfo> &sites,
                            const MetadataCache& cache,
                            const std::vector<std::string>& samples,
                            std::string &ynode);


class ResidualsFile {
private:
    std::string filename_;
    std::ofstream ofs_;

public:
    // constructor
    ResidualsFile(std::string filename) : filename_(filename) {}

    // destructor
    ~ResidualsFile();

    static Status Open(std::string filename, std::unique_ptr<ResidualsFile> &ans);

    // write a YAML record to the end of the file, as an element in a top-level sequence.
    Status write_record(std::string &rec);
};

} // namespace GLnexus
#endif
