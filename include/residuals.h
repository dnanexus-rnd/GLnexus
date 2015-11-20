#ifndef GLNEXUS_RESIDUALS_H
#define GLNEXUS_RESIDUALS_H

#include "data.h"
#include "types.h"
#include <fstream>
#include <memory>
#include "service_config.h"

namespace GLnexus {

// Write a record describing a loss into an output stream
Status write_residuals_record(std::ofstream &residuals_ofs,
                              const MetadataCache& cache, BCFData& data,
                              const unified_site& site,
                              const std::string& sampleset, const std::vector<std::string>& samples,
                              const bcf_hdr_t *gl_hdr, std::shared_ptr<bcf1_t>& gl_call);
} // namespace GLnexus
#endif
