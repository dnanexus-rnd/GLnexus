#include <assert.h>
#include "genotyper.h"

using namespace std;

namespace GLnexus {

// Placeholder using hard-called genotypes
Status translate_genotypes(const unified_site& site, const shared_ptr<const bcf_hdr_t>& dataset_header,
                           const shared_ptr<bcf1_t>& record, const map<int,int>& sample_mapping,
                           vector<int32_t>& genotypes, vector<bool>& genotyped) {
    assert(genotyped.size() > 0);
    assert(genotypes.size() == 2*genotyped.size());
    Status s;

    // map the BCF's alleles onto the unified alleles
    vector<int> allele_mapping;

    // reference allele maps if it contains the unified site
    range rng;
    S(range_of_bcf(dataset_header, record, rng));
    allele_mapping.push_back(rng.contains(site.pos) ? 0 : -1);

    // map alt alleles according to unification
    for (int i = 1; i < record->n_allele; i++) {
        auto p = site.unification.find(make_pair(rng.beg, string(record->d.allele[i])));
        allele_mapping.push_back(p != site.unification.end() ? p->second : -1);
    }

    // get the genotype calls
    int *gt = nullptr, gtsz = 0;
    int nGT = bcf_get_genotypes(dataset_header.get(), record.get(), &gt, &gtsz);
    assert(nGT == 2*bcf_hdr_nsamples(dataset_header.get()));

    // for each shared sample, record the genotype call.
    for (const auto& ij : sample_mapping) {
        assert(2*ij.first < nGT);
        assert(ij.second < genotyped.size());

        if (!genotyped[ij.second]) {
            // TODO: accept the call only if the allele is supported by at
            // least N reads -- including reference. The cutoff should be
            // configurable somehow.
            #define fill_allele(ofs)                                \
                if (gt[2*ij.first+ofs] != bcf_int32_vector_end &&   \
                    !bcf_gt_is_missing(gt[2*ij.first+ofs])) {       \
                    auto al = bcf_gt_allele(gt[2*ij.first+ofs]);    \
                    assert(al >= 0 && al < record->n_allele);       \
                    if (allele_mapping[al] >= 0) {                  \
                        genotypes[2*ij.second+ofs]                  \
                            = bcf_gt_unphased(allele_mapping[al]);  \
                    }                                               \
                }
            fill_allele(0)
            fill_allele(1)
            genotyped[ij.second] = true;
        } else {
            // subtlety: we already set the genotype for this
            // sample based on a previous BCF record. We now need
            // to set it back to missing, because we can't
            // accurately render the situation. TODO: support
            // symbolic non-ref allele
            genotypes[2*ij.second]
                = genotypes[2*ij.second+1]
                = bcf_gt_missing;
        }
    }

    if (gt) {
        free(gt);
    }

    return Status::OK();
}

Status genotype_site(const DataCache& data, const unified_site& site, const set<string>& samples,
					 const set<string>& datasets, const bcf_hdr_t* hdr, shared_ptr<bcf1_t>& ans) {
	Status s;
	
    // initialize a vector for the unified genotype calls (starting with everything missing)
    // TODO: ploidy
    vector<int32_t> genotypes(2*samples.size(), bcf_gt_missing);
    vector<bool> genotyped(samples.size(), false);

    // for each pertinent dataset
    for (const auto& dataset : datasets) {
        // load BCF records overlapping the site
        shared_ptr<const bcf_hdr_t> dataset_header;
        vector<shared_ptr<bcf1_t>> records;
        S(data.dataset_bcf(dataset, site.pos, dataset_header, records));

        // index the samples shared between the sample set and the BCFs
        // TODO: better algorithm, with caching (LRU given sampleset-dataset cross)
        map<int,int> sample_mapping;
        int bcf_nsamples = bcf_hdr_nsamples(dataset_header.get());
        for (int i = 0; i < bcf_nsamples; i++) {
            string sample_i(bcf_hdr_int2id(dataset_header.get(), BCF_DT_SAMPLE, i));
            int j = 0;
            for (const auto& sample_j : samples) {
                if (sample_i == sample_j) {
                    sample_mapping[i] = j;
                    break;
                }
                j++;
            }
        }

        // for each source BCF record
        for (const auto& record : records) {
            S(translate_genotypes(site, dataset_header, record, sample_mapping, genotypes, genotyped));
        }
    }

    // Create the destination BCF record for this site
    vector<const char*> c_alleles;
    for (const auto& allele : site.alleles) {
        c_alleles.push_back(allele.c_str());
    }

    ans = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);

    ans->rid = site.pos.rid;
    ans->pos = site.pos.beg;
    ans->rlen = site.pos.end - site.pos.beg;
    ans->qual = 0;
    if (bcf_update_alleles(hdr, ans.get(), c_alleles.data(), c_alleles.size()) != 0) {
        return Status::Failure("bcf_update_alleles");
    }
    if (bcf_update_genotypes(hdr, ans.get(), genotypes.data(), genotypes.size()) != 0) {
        return Status::Failure("bcf_update_genotypes");
    }

    return Status::OK();
}

}
