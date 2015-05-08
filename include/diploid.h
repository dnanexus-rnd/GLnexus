#ifndef GLNEXUS_DIPLOID_H
#define GLNEXUS_DIPLOID_H

#include <utility>

namespace diploid {

int genotypes(int n_allele);

std::pair<int,int> gt_alleles(int n_allele, int gt);
int alleles_gt(int n_allele, int a1, int a2);

namespace trio {
int mendelian_inconsistencies(int n_allele, int gt_p1, int gt_p2, int gt_ch);
}

}

#endif
