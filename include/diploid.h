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
// soft estimates of the allele copy numbers, assuming a uniform prior over
// genotypes.
// Otherwise, copy numbers are filled in based on the hard genotype calls.
Status estimate_allele_copy_number(const bcf_hdr_t* header, bcf1_t *record, std::vector<std::vector<float>>& ans);

namespace trio {
int mendelian_inconsistencies(int gt_p1, int gt_p2, int gt_ch);
}

}}


