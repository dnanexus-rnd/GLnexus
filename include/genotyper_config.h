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

    /// The symbolic allele used in INPUT gVCF reference confidence records
    std::string ref_symbolic_allele = "<NON_REF>";

    /// FORMAT field to consult for reference depth in gVCF reference records
    std::string ref_dp_format = "MIN_DP";

    /// Sometimes the input gVCF records contain specific allele calls that
    /// GLnexus is unable to render accurately in the output VCF file, e.g.
    /// rare alleles that are "pruned" during unification to ameliorate
    /// combinatorial explosion of the possible genotypes.
    /// By default, GLnexus emits no-calls (.) in such situations, which is
    /// the same thing that'd appear if the sample had no overlapping gVCF data
    /// at all. Sometimes it might be useful to distinguish these cases.
    /// If unrepr_symbolic_allele is set to a symbolic allele such as "?" or
    /// "NON_REF", then GLnexus will add this symbolic allele to the ALT of
    /// all output records, and calls of it indicate an unrepresentable allele
    /// call as described above.
    /// The symbolic allele must not contain angle brackets or whitespace.
    std::string unrepr_symbolic_allele;

    genotyper_config() = default;
};

}

#endif
