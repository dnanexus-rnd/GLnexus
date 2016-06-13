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

GLnexus::Status estimate_allele_copy_number(const bcf_hdr_t* header, bcf1_t *record, double bias00, vector<vector<double>>& ans) {
    ans.resize(record->n_sample);
    unsigned nGT = genotypes(record->n_allele);

    // buffers for genotype likelihoods
    htsvecbox<int32_t> igl;
    htsvecbox<float> fgl;
    vector<double> gl(record->n_sample*nGT);

    // try loading genotype likelihoods from PL
    if (bcf_get_format_int32(header, record, "PL", &igl.v, &igl.capacity) == record->n_sample*nGT) {
        for (unsigned ik = 0; ik < record->n_sample*nGT; ik++) {
            assert(igl[ik] >= 0.0);
            gl[ik] = exp10(double(igl[ik])/(-10.0));
        }
/*
    // couldn't load PL; try GL
    } else if (bcf_get_format_float(header, record, "GL", &fgl.v, &fgl.capacity) == record->n_sample*nGT) {
        for (unsigned ik = 0; ik < record->n_sample*nGT; ik++) {
            assert(fgl[ik] <= 0.0);
            gl[ik] = exp10(gl[ik]));
        } */
    } else {
        gl.clear();
    }

    if (!gl.empty()) {
        for (unsigned i = 0; i < record->n_sample; i++) {
            double *gl_i = &(gl[i*nGT]);
            vector<double>& ans_i = ans[i];
            ans_i.resize(record->n_allele);
            memset(ans_i.data(), 0, record->n_allele*sizeof(double));

            // based on genotype likelihoods and specified bias,
            // calculate expected # of copies of each allele
            double total_likelihood = 0.0;
            for (unsigned k = 0; k < nGT; k++) {
                double gl_ik = gl_i[k] * (k == 0 ? bias00 : 1.0);
                total_likelihood += gl_ik;
                auto p = gt_alleles(k);
                ans_i[p.first] += gl_ik;
                ans_i[p.second] += gl_ik;
            }
            for (unsigned j = 0; j < record->n_allele; j++) {
                ans_i[j] /= total_likelihood;
            }
        }
    } else {
        // couldn't find any genotype likelihoods; set copy number based on
        // the hard genotype calls
        htsvecbox<int> gt;
        if (bcf_get_genotypes(header, record, &gt.v, &gt.capacity) != 2*record->n_sample) {
            return Status::Invalid("diploid::estimate_copy_number: couldn't load PL, GL, nor genotypes", range(record).str());
        }

        for (unsigned i = 0; i < record->n_sample; i++) {
            vector<double>& ans_i = ans[i];
            ans_i.resize(record->n_allele);
            memset(ans_i.data(), 0, record->n_allele*sizeof(double));
            for (unsigned ofs : {0, 1}) {
                if (!bcf_gt_is_missing(gt[i*2+ofs])) {
                    unsigned j = bcf_gt_allele(gt[i*2+ofs]);
                    assert(j < record->n_allele);
                    ans_i[j] += 1.0;
                }
            }
        }
    }

    return Status::OK();
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
