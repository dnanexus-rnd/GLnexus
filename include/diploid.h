#pragma once

#include <utility>
#include <vector>
#include "types.h"

namespace GLnexus {
namespace diploid {

unsigned genotypes(unsigned n_allele);

std::pair<unsigned,unsigned> gt_alleles(unsigned gt);
unsigned alleles_gt(unsigned a1, unsigned a2);

// Estimate allele copy number for each sample from a BCF record.
// Produces a record->n_sample*record->n_allele matrix with an estimate, for
// each sample, of the copy number of each allele.
// If PL or GL genotype likelihoods are available, these are used to derive
// soft estimates of the allele copy numbers. A prior distribution on the genotypes
// is applied in deriving these estimates, flat across all genotypes except
// homozygous ref (0/0), which is multiplied by bias00. That is, the prior over
// genotypes [0/0, 0/1, 1/1, ...] is proportional to [bias00, 1, 1, ..].
// Otherwise, copy numbers are filled in based on the hard genotype calls.
Status estimate_allele_copy_number(const bcf_hdr_t* header, bcf1_t *record, double bias00, std::vector<std::vector<double>>& ans);
Status bcf_alleles_maxAQ(const bcf_hdr_t* hdr, bcf1_t* record, const std::vector<unsigned>& samples, std::vector<int>& ans);

namespace trio {
int mendelian_inconsistencies(int gt_p1, int gt_p2, int gt_ch);
}

}}


