#include <utility>
#include <math.h>
#include <assert.h>
#include <sstream>
#include "diploid.h"
using namespace std;

namespace GLnexus {
namespace diploid {

// n_gt = (nA+1) choose 2 = (nA+1)!/2/(nA-1)! = (nA+1)(nA)/2
unsigned genotypes(unsigned n_allele) {
    return (n_allele+1)*n_allele/2;
}

unsigned alleles_gt(unsigned i, unsigned j) {
    if (j < i) {
        swap(i, j);
    }
    return (j*(j+1)/2)+i;
}

// given a genotype index, return a pair with the indices of the constituent n_allele
// TODO: memoize for small values of n_allele...
pair<unsigned,unsigned> gt_alleles(unsigned gt) {
    /*
      0 1 2
    0 0 1 3
    1   2 4
    2     5

      0 1 2 3
    0 0 1 3 6
    1   2 4 7
    2     5 8
    3       9
    */

    unsigned j = (unsigned)((sqrt(8*gt+1)-1.0)/2.0);
    assert(gt >= j*(j+1)/2);
    unsigned i = gt - j*(j+1)/2;
    return make_pair(i,j);
}

Status bcf_zygosity_by_GQ(const bcf_hdr_t* header, bcf1_t* record, const std::vector<unsigned>& samples,
                          vector<zygosity_by_GQ>& ans) {
    htsvecbox<int> gt;
    int nGT = bcf_get_genotypes(header, record, &gt.v, &gt.capacity);
    if (record->n_sample == 1 && nGT == 1 && !gt.empty()) {
        // special case for Strelka2 and other callers which emit some gVCF
        // records with GT=. or GT=0 or GT=1: rewrite these to look like ./.
        // and ./0 and ./1 as far as our genotyper is concerned.
        gt.v = (int*) realloc(gt.v, 2*sizeof(int));
        gt.capacity = 2;
        swap(gt[0], gt[1]);
        gt[0] = bcf_gt_missing;
        assert(bcf_gt_is_missing(gt[0]));
        assert(bcf_gt_is_missing(gt[1]) || bcf_gt_allele(gt[1]) == 0);
        nGT = 2;
    }
    if (gt.empty() || nGT != 2*record->n_sample) return Status::Failure("bcf_get_genotypes");

    htsvecbox<int32_t> gq;
    int nGQ = bcf_get_format_int32(header, record, "GQ", &gq.v, &gq.capacity);
    if (nGQ != record->n_sample) {
        gq.clear();
    }

    ans.resize(record->n_allele);
    for (auto& z : ans) {
        z.clear();
    }

    for (auto sample : samples) {
        auto sample_gq = gq.empty() ? 0 : gq[sample];
        auto pres1 = !bcf_gt_is_missing(gt[sample*2]), pres2 = !bcf_gt_is_missing(gt[sample*2+1]);
        auto al1 = pres1 ? bcf_gt_allele(gt[sample*2]) : -1;
        auto al2 = pres2 ? bcf_gt_allele(gt[sample*2+1]) : -1;

        if (pres1 && pres2 && al1 == al2) {
            assert(al1 >= 0 && al1 < record->n_allele);
            ans[al1].add(2,sample_gq);
            assert(ans[al1].copy_number() > 0);
        } else {
            if (pres1) {
                assert(al1 >= 0 && al1 < record->n_allele);
                ans[al1].add(1,sample_gq);
                assert(ans[al1].copy_number() > 0);
            }
            if (pres2) {
                assert(al2 >= 0 && al2 < record->n_allele);
                ans[al2].add(1,sample_gq);
                assert(ans[al2].copy_number() > 0);
            }
        }
    }

    return Status::OK();
}

const double LOG0 = log(0.0);
const double LOG10_E = log10(exp(1.0));
GLnexus::Status bcf_get_genotype_log_likelihoods(const bcf_hdr_t* header, bcf1_t *record, vector<double>& gll) {
    unsigned nGT = genotypes(record->n_allele);
    gll.resize(record->n_sample*nGT);

    // try loading genotype likelihoods from PL
    htsvecbox<int32_t> igl;
    if (bcf_get_format_int32(header, record, "PL", &igl.v, &igl.capacity) == record->n_sample*nGT) {
        for (unsigned ik = 0; ik < record->n_sample*nGT; ik++) {
            auto x = igl[ik];
            if (x == bcf_int32_missing || x == bcf_int32_vector_end) {
                gll[ik] = LOG0;
            } else if (x >= 0) {
                gll[ik] = double(igl[ik])/(-10.0*LOG10_E);
            } else {
                return Status::Invalid("bcf_get_genotype_log_likelihoods: negative PL entry");
            }
        }
        return Status::OK();
    }
/*
    // couldn't load PL; try GL
    htsvecbox<float> fgl;
    if (bcf_get_format_float(header, record, "GL", &fgl.v, &fgl.capacity) == record->n_sample*nGT) {
        for (unsigned ik = 0; ik < record->n_sample*nGT; ik++) {
            assert(fgl[ik] <= 0.0);
            gl[ik] = double(gl[ik])/LOG10_E;
        }
        return Status::OK();
    }
*/
    return Status::NotFound();
}

// for each allele, find the max AQs across the given sample indices, given valid genotype log-likelihoods
GLnexus::Status alleles_topAQ(unsigned n_allele, unsigned n_sample, const vector<unsigned>& samples,
                              const vector<double>& gll, vector<top_AQ>& ans) {
    unsigned nGT = genotypes(n_allele);
    if (gll.size() != n_sample*nGT) return Status::Invalid("alleles_topAQ");
    vector<vector<int>> obs(n_allele);

    for (unsigned i : samples) {
        assert(i < n_sample);
        const double *gll_i = gll.data() + i*nGT;
        vector<double> maxLL_with(n_allele, LOG0);    // max likelihood of a genotype carrying each allele
        vector<double> maxLL_without(n_allele, LOG0); // max likelihood of a genotype NOT carrying each allele

        for (unsigned k = 0; k < nGT; k++) { // for each genotype
            auto p = gt_alleles(k);          // indices of the two alleles in this genotype
            for (unsigned al = 0; al < n_allele; al++) {
                if (al == p.first || al == p.second) {
                    // genotype contains this allele -- update maxL_with
                    maxLL_with[al] = std::max(maxLL_with[al], gll_i[k]);
                } else {
                    // genotype doesn't contain this allele -- update maxL_without
                    maxLL_without[al] = std::max(maxLL_without[al], gll_i[k]);
                }
            }
        }
        for (unsigned al = 0; al < n_allele; al++) {
            if (maxLL_with[al] == LOG0) {
                obs[al].push_back(0);
            } else if (maxLL_without[al] == LOG0) {
                obs[al].push_back(MAX_AQ);
            } else {
                // phred scale likelihood ratio
                double AQLLR = std::max(0.0, maxLL_with[al] - maxLL_without[al]);
                obs[al].push_back(std::min(MAX_AQ,(int)round(10.0*AQLLR/log(10.0))));
            }
        }
    }

    ans.resize(n_allele);
    for (unsigned i = 0; i < n_allele; i++) {
        ans[i].clear();
        ans[i] += obs[i];
    }

    return Status::OK();
}

// for each allele, find the top AQs across the given sample indices, from the BCF record
// if no genotype likelihoods can be found in the record, then return zeroes
GLnexus::Status bcf_alleles_topAQ(const bcf_hdr_t* hdr, bcf1_t* record, const vector<unsigned>& samples,
                                  vector<top_AQ>& ans) {
    vector<double> gll;
    GLnexus::Status s = bcf_get_genotype_log_likelihoods(hdr, record, gll);

    if (s.ok()) {
        return alleles_topAQ(record->n_allele, record->n_sample, samples, gll, ans);
    } else if (s == GLnexus::StatusCode::NOT_FOUND) {
        ans.resize(record->n_allele);
        for (auto& v : ans) {
            v.clear();
        }
        return Status::OK();
    }

    return s;
}

namespace trio {
// Given genotype indices of two parents and one child, return:
//   0 if one of the child alleles is present in one parent and the other child
//   allele is present in the other parent
//   1 if one child allele is found in a parent, but the other child allele is
//   not found in the other parent
//   2 if neither child allele is found in the parents
// Note, this is not simply number of alleles not found in either parent, e.g.
// p1=0/0, p2=1/1, ch=1/1 => 1
int mendelian_inconsistencies(int gt_p1, int gt_p2, int gt_ch) {
    auto a_p1 = gt_alleles(gt_p1);
    auto a_p2 = gt_alleles(gt_p2);
    auto a_ch = gt_alleles(gt_ch);
    auto a1 = a_ch.first, a2 = a_ch.second;

    if (a1 == a_p1.first || a1 == a_p1.second) {
        return (a2 == a_p2.first || a2 == a_p2.second) ? 0 : 1;
    } else if (a1 == a_p2.first || a1 == a_p2.second) {
        return (a2 == a_p1.first || a2 == a_p1.second) ? 0 : 1;
    } else {
        return (a2 == a_p1.first || a2 == a_p1.second || a2 == a_p2.first || a2 == a_p2.second) ? 1 : 2;
    }
}
}}}
