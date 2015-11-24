#ifndef GLNEXUS_RESIDUALS_H
#define GLNEXUS_RESIDUALS_H

#include <mutex>
#include <fstream>
#include <memory>
#include "data.h"
#include "types.h"
#include "service_config.h"

namespace GLnexus {

class Residuals {
private:
    std::string filename_; // file where to store the residual calls
    std::mutex mutex_;
    std::ofstream ofs_;
    const MetadataCache& cache_;
    BCFData& data_;
    const std::string& sampleset_;
    const std::vector<std::string>& samples_;

    // non lock safe version
    Status write_record_(const unified_site& site,
                         const bcf_hdr_t *gl_hdr,
                         bcf1_t *gl_call);

public:
    // constructor
    Residuals(std::string filename,
              const MetadataCache& cache, BCFData& data,
              const std::string& sampleset, const std::vector<std::string>& samples) :
        filename_(filename), cache_(cache), data_(data), sampleset_(sampleset), samples_(samples)
        {}

    // destructor
    ~Residuals();

    // Create a Residuals object
    static Status Open(std::string filename,
                       const MetadataCache& cache, BCFData& data,
                       const std::string& sampleset, const std::vector<std::string>& samples,
                       std::unique_ptr<Residuals> &ans);

    // Write a record describing a loss
    //
    // note: this function is called from multiple threads.
    Status write_record(const unified_site& site,
                        const bcf_hdr_t *gl_hdr,
                        bcf1_t *gl_call);
};

} // namespace GLnexus
#endif
