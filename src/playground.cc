#include <iostream>
#include <vcf.h>
#include <vector>
#include <limits>
#include <memory>
#include <math.h>
#include <assert.h>
#include "diploid.h"
using namespace std;

extern "C" int bcf_hdr_register_hrec(bcf_hdr_t *hdr, bcf_hrec_t *hrec);

// sugar for declaring a unique_ptr with a custom deleter function
#define UPD(T) std::unique_ptr<T, void(*)(T*)>

int main() {
    auto vcf = bcf_open("test/data/trio_denovo.vcf", "r");

    UPD(bcf_hdr_t) hdr(bcf_hdr_read(vcf), &bcf_hdr_destroy);

    bcf_hrec_t chr21;
    chr21.type = BCF_HL_CTG;
    chr21.key = "contig";
    chr21.value = nullptr;
    chr21.nkeys = 2;
    unique_ptr<char*[]> chr21keys(new char*[2]);
    chr21keys[0] = "ID";
    chr21keys[1] = "length";
    chr21.keys = chr21keys.get();
    unique_ptr<char*[]> chr21vals(new char*[2]);
    chr21vals[0] = "21";
    chr21vals[1] = "10000000";
    chr21.vals = chr21vals.get();

    bcf_hdr_register_hrec(hdr.get(), &chr21);
    bcf_hdr_sync(hdr.get());


    // TODO augment header with hs37d5 sequences in the BCF_DT_CTG
    // -D FILE   sequence dictionary for VCF->BCF conversion [null]
/*
    for (int i = 0; i < hdr->nhrec; i++) {
        auto key = hdr->hrec[i]->key;
        auto value = hdr->hrec[i]->value;

        if (value) {
            cout << key << "=" << value << endl;
        }
    }

    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
        cout << hdr->samples[i] << " ";
    }
    cout << endl;
*/
    auto vt = bcf_init();

    int *gt = nullptr, gtsz = 0, glsz = 0, nGT, nGL;
    float *gl = nullptr;

    auto stdout = vcf_open("-", "w");

    bcf_hdr_write(stdout, hdr.get());

    while (bcf_read(vcf, hdr.get(), vt) == 0) {
        vt->errcode &= ~BCF_ERR_CTG_UNDEF;
        bcf_unpack(vt, BCF_UN_ALL);

        bcf_write(stdout, hdr.get(), vt);

        nGL = bcf_get_format_float(hdr.get(), vt, "GL", &gl, &glsz);

        int n_gt = GLnexus::diploid::genotypes(vt->n_allele);
        float max_post = -std::numeric_limits<float>::infinity();
        int fa_max, mo_max, pb_max, sb_max;

        // TODO can this be expressed in a clever vectorized form?
        for (int fa = 0; fa < n_gt; fa++) {
            for (int mo = 0; mo < n_gt; mo++) {
                for (int pb = 0; pb < n_gt; pb++) {
                    for (int sb = 0; sb < n_gt; sb++) {
                        int denovo = GLnexus::diploid::trio::mendelian_inconsistencies(fa, mo, pb) +
                                     GLnexus::diploid::trio::mendelian_inconsistencies(fa, mo, sb);

                        float post =  gl[fa] + gl[3+mo] + gl[6+pb] + gl[9+sb] - 6*denovo;
                        if (post > max_post) {
                            max_post = post;
                            fa_max = fa; mo_max = mo;
                            pb_max = pb; sb_max = sb;
                        }
                    }
                }
            }
        }

        nGT = bcf_get_genotypes(hdr.get(), vt, &gt, &gtsz);

        for (int i = 0; i < nGT; i++) {
            gt[i] = bcf_gt_unphased(0);
        }
        #define SET_GT(map_gt, ofs) \
                    if (map_gt>1) gt[ofs] = bcf_gt_unphased(1); \
                    if (map_gt) gt[ofs+1] = bcf_gt_unphased(1);
        SET_GT(fa_max, 0)
        SET_GT(mo_max, 2)
        SET_GT(pb_max, 4)
        SET_GT(sb_max, 6)
        // manipulate alleles, likelihoods, ...
        bcf_update_genotypes(hdr.get(), vt, gt, nGT);

        bcf_write(stdout, hdr.get(), vt);
    }

    vcf_close(stdout);

    if (gt) free(gt);
    if (gl) free(gl);

    bcf_destroy(vt);
    bcf_close(vcf);

    return 0;
}
