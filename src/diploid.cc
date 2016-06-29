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
    if (gt.empty() || nGT != 2*record->n_sample) return Status::Failure("bcf_get_genotypes");

    htsvecbox<int32_t> gq;
    int nGQ = bcf_get_format_int32(header, record, "GQ", &gq.v, &gq.capacity);
    if (nGQ != record->n_sample) {
        gq.clear();
    }

    ans.clear();
    ans.resize(record->n_allele);

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

GLnexus::Status bcf_get_genotype_likelihoods(const bcf_hdr_t* header, bcf1_t *record, vector<double>& gl) {
    unsigned nGT = genotypes(record->n_allele);
    gl.resize(record->n_sample*nGT);

     // try loading genotype likelihoods from PL
    htsvecbox<int32_t> igl;
    if (bcf_get_format_int32(header, record, "PL", &igl.v, &igl.capacity) == record->n_sample*nGT) {
        for (unsigned ik = 0; ik < record->n_sample*nGT; ik++) {
            assert(igl[ik] >= 0.0);
            gl[ik] = exp10(double(igl[ik])/(-10.0));
        }
        return Status::OK();
    }
/*
    // couldn't load PL; try GL
    htsvecbox<float> fgl;
    if (bcf_get_format_float(header, record, "GL", &fgl.v, &fgl.capacity) == record->n_sample*nGT) {
        for (unsigned ik = 0; ik < record->n_sample*nGT; ik++) {
            assert(fgl[ik] <= 0.0);
            gl[ik] = exp10(gl[ik]);
        }
        return Status::OK();
    }
*/
    return Status::NotFound();
}

// for each allele, find the max AQ across the given sample indices, given valid genotype likelihoods
GLnexus::Status alleles_maxAQ(unsigned n_allele, unsigned n_sample, const vector<unsigned>& samples,
                              const vector<double>& gl, vector<int>& ans) {
    unsigned nGT = genotypes(n_allele);
    if (gl.size() != n_sample*nGT) return Status::Invalid("alleles_maxAQ");
    ans.assign(n_allele,0);

    for (unsigned i : samples) {
        assert(i < n_sample);
        const double *gl_i = gl.data() + i*nGT;
        vector<double> maxL_with(n_allele, 0.0);    // max likelihood of a genotype carrying each allele
        vector<double> maxL_without(n_allele, 0.0); // max likelihood of a genotype NOT carrying each allele

        for (unsigned k = 0; k < nGT; k++) { // for each genotype
            auto p = gt_alleles(k);          // indices of the two alleles in this genotype
            for (unsigned al = 0; al < n_allele; al++) {
                if (al == p.first || al == p.second) {
                    // genotype contains this allele -- update maxL_with
                    maxL_with[al] = std::max(maxL_with[al], gl_i[k]);
                } else {
                    // genotype doesn't contain this allele -- update maxL_without
                    maxL_without[al] = std::max(maxL_without[al], gl_i[k]);
                }
            }
        }
        for (unsigned al = 0; al < n_allele; al++) {
            // phred scale likelihood ratio
            double AQLR = std::max(1.0, maxL_with[al]/maxL_without[al]);
            ans[al] = std::max(ans[al], (int)round(10.0*log(AQLR)/log(10.0)));
        }
    }

    return Status::OK();
}

// for each allele, find the max AQ across the given sample indices, from the BCF record
// if no genotype likelihoods can be found in the record, then return [0, 0, ...]
GLnexus::Status bcf_alleles_maxAQ(const bcf_hdr_t* hdr, bcf1_t* record, const vector<unsigned>& samples, vector<int>& ans) {
    ans.assign(record->n_allele,0);
    vector<double> gl;
    GLnexus::Status s = bcf_get_genotype_likelihoods(hdr, record, gl);

    if (s.ok()) {
        return alleles_maxAQ(record->n_allele, record->n_sample, samples, gl, ans);
    } else if (s == GLnexus::StatusCode::NOT_FOUND) {
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
