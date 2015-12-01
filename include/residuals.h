#ifndef GLNEXUS_RESIDUALS_H
#define GLNEXUS_RESIDUALS_H

#include <fstream>
#include <memory>
#include "data.h"
#include "types.h"
#include "service_config.h"

namespace GLnexus {

class Residuals {
private:
    const MetadataCache& cache_;
    BCFData& data_;
    const std::string& sampleset_;
    const std::vector<std::string>& samples_;

public:
    // constructor
    Residuals(const MetadataCache& cache, BCFData& data,
              const std::string& sampleset, const std::vector<std::string>& samples) :
        cache_(cache), data_(data), sampleset_(sampleset), samples_(samples)
        {}

    // destructor
    ~Residuals();

    // Create a Residuals object
    static Status Open(const MetadataCache& cache, BCFData& data,
                       const std::string& sampleset, const std::vector<std::string>& samples,
                       std::unique_ptr<Residuals> &ans);

    // Create a YAML node describing a loss. The node
    // is formatted as a string for simplicity.
    Status gen_record(const unified_site& site,
                      const bcf_hdr_t *gl_hdr,
                      const bcf1_t *gl_call,
                      std::string &ynode);
};

} // namespace GLnexus
#endif
