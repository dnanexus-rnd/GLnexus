#ifndef GLNEXUS_GENOTYPER_CONFIG_H
#define GLNEXUS_GENOTYPER_CONFIG_H

#include <types.h>
#include <string>

namespace GLnexus {

struct genotyper_config {
    /// Require any allele call to be supported by at least this depth
    size_t required_dp = 0;

    /// FORMAT field to consult for per-allele depth in VCF records
    std::string allele_dp_format = "AD";

    /// The symbolic allele used in gVCF reference confidence models
    std::string ref_symbolic_allele = "<NON_REF>";

    /// FORMAT field to consult for reference depth in gVCF reference records
    std::string ref_dp_format = "MIN_DP";

    genotyper_config() = default;
};

}

#endif
