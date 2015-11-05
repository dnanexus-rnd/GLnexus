#pragma once

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

    /// Output format (default = bcf), choices = "BCF", "VCF"
    std::string output_format = "BCF";

    genotyper_config() = default;

    genotyper_config(std::string _output_format) : output_format(_output_format) {}
};

}

