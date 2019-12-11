#include <vcf.h>

#include <string>
#include <vector>

#include "data.h"
#include "types.h"

namespace GLnexus {

// test whether a gVCF file is compatible for deposition into the database.
bool gvcf_compatible(const MetadataCache& metadata, const bcf_hdr_t* hdr);

// hard-coded xAtlas ingestion exceptions
// filter VRFromDeletion: accessory information
// format VR+RR >= 65536: unreliable QC values
bool xAtlas_ingestion_exceptions(const bcf_hdr_t* hdr, bcf1_t* bcf);

// Sanity-check an individual bcf1_t record before ingestion.
Status validate_bcf(const std::vector<std::pair<std::string, size_t>>& contigs,
                    const std::string& filename, const bcf_hdr_t* hdr,
                    bcf1_t* bcf, int prev_rid, int prev_pos,
                    bool& skip_ingestion);

// Verify that a VCF file is well formed.
// AND, fill in the [samples_out]
Status vcf_validate_basic_facts(MetadataCache& metadata,
                                const std::string& dataset,
                                const std::string& filename, bcf_hdr_t* hdr,
                                vcfFile* vcf,
                                std::set<std::string>& samples_out);

}  // namespace GLnexus
