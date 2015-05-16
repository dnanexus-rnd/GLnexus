#ifndef GLNEXUS_GENOTYPER_H
#define GLNEXUS_GENOTYPER_H

#include "data.h"
#include <memory>

namespace GLnexus {

Status genotype_site(const DataCache& data, const unified_site& site,
					 const std::set<std::string>& samples, const std::set<std::string>& datasets,
					 const bcf_hdr_t* hdr, std::shared_ptr<bcf1_t>& ans);

}

#endif
