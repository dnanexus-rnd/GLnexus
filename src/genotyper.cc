#include <assert.h>
#include "genotyper.h"
#include <iostream>

using namespace std;

// Here we implement a placeholder algorithm for genotyping individual samples
// at unified sites. It's merely capable of substituting in hard genotype
// calls for exactly matching alleles from our gVCF input data (with some
// simple filters on depth of coverage). It does not handle genotype
// likelihoods or reference base padding. Also it assumes diploid.

namespace GLnexus {

bool is_gvcf_ref_record(const genotyper_config& cfg, const bcf1_t* record) {
    return record->n_allele == 2 && string(record->d.allele[1]) == cfg.ref_symbolic_allele;
}

// Helper class for keeping track of the per-allele depth of coverage info in
// a bcf1_t record. There are a couple different cases to handle, depending on
// whether we're looking at a gVCF reference confidence record or a "regular"
// VCF record.
class AlleleDepthHelper {
    const genotyper_config& cfg_;
    size_t n_sample_, n_allele_;
    bool is_g_;
    shared_ptr<int32_t> v_;

    AlleleDepthHelper(const genotyper_config& cfg,
                      size_t n_sample, size_t n_allele, bool is_g,
                      const shared_ptr<int32_t>& v)
        : cfg_(cfg), n_sample_(n_sample), n_allele_(n_allele), is_g_(is_g), v_(v)
        {}

public:
    // TODO refactor to avoid heap allocation...
    static Status Open(const genotyper_config& cfg, const string& dataset,
                       const bcf_hdr_t* dataset_header, bcf1_t* record,
                       unique_ptr<AlleleDepthHelper>& ans) {
        if (cfg.required_dp) {
            int32_t *v = nullptr;
            int vsz = 0;
            bool is_g = false;

            // is this a gVCF reference confidence record?
            if (is_gvcf_ref_record(cfg, record)) {
                is_g = true;
                // if so, look for the MIN_DP FORMAT field (or equivalent)
                int nv = bcf_get_format_int32(dataset_header, record, cfg.ref_dp_format.c_str(),
                                              &v, &vsz);

                if (nv != record->n_sample) {
                    if (v) free(v);
                    ostringstream errmsg;
                    errmsg << dataset << " " << range(record).str() << " (" << cfg.ref_dp_format << ")";
                    return Status::Invalid("genotyper: gVCF reference depth FORMAT field is either missing or has the wrong type", errmsg.str());
                }
            } else {
                // this is a regular VCF record, so look for the AD FORMAT field (or equivalent)
                int nv = bcf_get_format_int32(dataset_header, record, cfg.allele_dp_format.c_str(),
                                              &v, &vsz);

                if (nv != record->n_sample * record->n_allele) {
                    if (v) free(v);
                    ostringstream errmsg;
                    errmsg << dataset << " " << range(record).str() << " (" << cfg.allele_dp_format << ")";
                    return Status::Invalid("genotyper: VCF allele depth FORMAT field is either missing or has the wrong type", errmsg.str());
                }
            }

            ans.reset(new AlleleDepthHelper(cfg, size_t(record->n_sample), size_t(record->n_allele), is_g,
                                            shared_ptr<int32_t>(v, [](int32_t* v) { free(v); })));
        } else {
            ans.reset(new AlleleDepthHelper(cfg, 0, 0, false, nullptr));
        }
        return Status::OK();
    }

    // Is there sufficient coverage for sample i, allele j?
    bool sufficient(size_t sample, size_t allele) {
        if (!cfg_.required_dp) {
            return true;
        }
        assert(sample < n_sample_);
        assert(allele < n_allele_);
        if (is_g_) {
            // the MIN_DP array has just one integer per sample
            if (allele == 0) {
                return v_.get()[sample] >= cfg_.required_dp;
            } else {
                return false;
            }
        } else {
            // the AD array has one integer per allele per sample
            return v_.get()[sample*n_allele_+allele] >= cfg_.required_dp;
        }
    }
};

// Translate the hard-called genotypes from the bcf1_t into our genotype
// vector, based on a mapping from the bcf1_t sample indices into indices of
// the genotype vector.
Status translate_genotypes(const genotyper_config& cfg, const unified_site& site,
                           const string& dataset, const bcf_hdr_t* dataset_header,
                           bcf1_t* record, const map<int,int>& sample_mapping,
                           vector<int32_t>& genotypes, vector<bool>& genotyped) {
    assert(genotyped.size() > 0);
    assert(genotypes.size() == 2*genotyped.size());
    Status s;

    // map the BCF's alleles onto the unified alleles
    vector<int> allele_mapping;

    // reference allele maps if it contains the unified site
    range rng(record);
    allele_mapping.push_back(rng.contains(site.pos) ? 0 : -1);

    // map the bcf1_t alt alleles according to unification
    for (int i = 1; i < record->n_allele; i++) {
        auto p = site.unification.find(make_pair(rng.beg, string(record->d.allele[i])));
        allele_mapping.push_back(p != site.unification.end() ? p->second : -1);
    }

    // get the genotype calls
    int *gt = nullptr, gtsz = 0;
    int nGT = bcf_get_genotypes(dataset_header, record, &gt, &gtsz);
    assert(nGT == 2*bcf_hdr_nsamples(dataset_header));
    assert(record->n_sample == bcf_hdr_nsamples(dataset_header));

    // and the depth of coverage info
    unique_ptr<AlleleDepthHelper> depth;
    S(AlleleDepthHelper::Open(cfg, dataset, dataset_header, record, depth));

    // for each shared sample, record the genotype call.
    for (const auto& ij : sample_mapping) {
        assert(2*ij.first < nGT);
        assert(ij.second < genotyped.size());

        if (!genotyped[ij.second]) {
            #define fill_allele(ofs)                                \
                if (gt[2*ij.first+ofs] != bcf_int32_vector_end &&   \
                    !bcf_gt_is_missing(gt[2*ij.first+ofs])) {       \
                    auto al = bcf_gt_allele(gt[2*ij.first+ofs]);    \
                    assert(al >= 0 && al < record->n_allele);       \
                    if (allele_mapping[al] >= 0 && depth->sufficient(ij.second, al)) { \
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
            // accurately render the situation.

            int lost_index = 2*ij.second;
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

Status update_orig_loss(const bcf_hdr_t* dataset_header, bcf1_t* record, map<int,int>& sample_mapping, vector<loss_stats>& losses) {
    range rng(record);
    Status s;

    int *gt = nullptr, gtsz = 0;
    int nGT = bcf_get_genotypes(dataset_header, record, &gt, &gtsz);
    int n_bcf_samples = bcf_hdr_nsamples(dataset_header);
    assert(nGT == 2*n_bcf_samples);
    assert(record->n_sample == n_bcf_samples);

    for (int i = 0; i < n_bcf_samples; i++) {
        int sample_ind = sample_mapping[i];
        auto& loss = losses[sample_ind];

        // No calls on this sample, move to next sample 
        if(bcf_gt_is_missing(gt[i*2]) && bcf_gt_is_missing(gt[i*2 + 1]))
            continue;

        // Check if the call range already exists in the list of calls recorded
        // TODO: handling of overlapping calls
        if (loss.orig_calls.find(rng) == loss.orig_calls.end())
            loss.orig_calls.insert(rng);
        else
            continue;

        if (!bcf_gt_is_missing(gt[i*2])){
            // cout << "Adding to " << loss.sample_name << " disp 0; range " << rng.str() << " " << rng.size() << endl;
            loss.orig_bp += rng.size();
        }
        if (!bcf_gt_is_missing(gt[i*2 + 1])){
            // cout << "Adding to " << loss.sample_name << " disp 1; range " << rng.str() << " " << rng.size()  << endl;
            loss.orig_bp += rng.size();
        }
    }

    return Status::OK();
}

Status update_called_loss(const bcf_hdr_t* hdr, bcf1_t* record, vector<loss_stats>& losses) {
    range rng(record);
    Status s;

    int *gt = nullptr, gtsz = 0;
    int nGT = bcf_get_genotypes(hdr, record, &gt, &gtsz);
    int n_bcf_samples = record->n_sample;
    assert(n_bcf_samples == losses.size());

    for (int i = 0; i < n_bcf_samples; i++) {
        auto& loss = losses[i];
        
        // Assumes no overlapping calls/ calls in the same range in the output joint calls
        if(!bcf_gt_is_missing(gt[i*2]) || !bcf_gt_is_missing(gt[i*2 + 1]))
            loss.joint_calls.insert(rng);

        if (!bcf_gt_is_missing(gt[i*2]))
            loss.joint_bp += rng.size();
        if (!bcf_gt_is_missing(gt[i*2 + 1]))
            loss.joint_bp += rng.size();
    }

    return Status::OK();
}   

Status genotype_site(const genotyper_config& cfg, const DataCache& data, const unified_site& site,
                     const set<string>& samples, const set<string>& datasets,
                     const bcf_hdr_t* hdr, shared_ptr<bcf1_t>& ans, vector<loss_stats>& losses) {
	Status s;
	
    // Initialize a vector for the unified genotype calls for each sample,
    // starting with everything missing. We'll then loop through BCF records
    // overlapping this site and fill in the genotypes as we encounter them.
    vector<int32_t> genotypes(2*samples.size(), bcf_gt_missing);
    // Also remember which samples we've already seen a genotype call for, in
    // case we encounter multiple BCF records from the sample
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
            S(translate_genotypes(cfg, site, dataset, dataset_header.get(), record.get(),
                                  sample_mapping, genotypes, genotyped));

            S(update_orig_loss(dataset_header.get(), record.get(), sample_mapping, losses));
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

    S(update_called_loss(hdr, ans.get(), losses));

    return Status::OK();
}

}
