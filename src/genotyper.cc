#include <assert.h>
#include <algorithm>
#include "genotyper.h"

using namespace std;

// Here we implement a v1 early algorithm for genotyping individual samples
// at unified sites. It's currently capable of substituting in hard genotype
// calls for exactly matching alleles from our gVCF input data, with reference
// padding using neighboring gvcf reference records (where necessary). It does
// not handle joining of multiple variant records within a single site.
// Basic filtering based on coverage is supported.
// It does not handle genotype likelihoods, or carry over other fields besides
// GT (and RNC for accountability). Also it assumes diploid.
namespace GLnexus {

/// Determine whether the given record is a gVCF reference confidence record
/// (or else a "normal" record with at least one specific ALT allele)
static bool is_gvcf_ref_record(const genotyper_config& cfg, const bcf1_t* record) {
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

    // Whether there is sufficient ref coverage, expect depth
    // passed in to be < 0 if the representative record had no reference
    // region
    bool sufficient_ref(size_t depth) {
        // Return true if no required_dp, or record has no ref region
        if (!cfg_.required_dp || depth < 0) {
            return true;
        }
        return (depth >= cfg_.required_dp);
    }

    // Returns the depth of the gvcf_ref coverage for sample i
    // returns -1 if record is not a gvcf confidence record
    int get_gvcf_depth(size_t sample) {
        if (!is_g_) {
            return -1;
        } else {
            return v_.get()[sample];
        }
    }
};


// Called for each bcf record that is associated with the unified_site
// examined to update "denominator" of total calls/bp covered in orig
// dataset
Status LossTracker::add_call_for_site(const range call, int n_calls, bool is_gvcf) noexcept {
    if (is_finalized)
        return Status::Invalid("calling add_call_for_site for a finalized LossTracker");

    auto rng_within_site_p = call.intersect(rng);

    if (rng_within_site_p) {
        range rng_within_site = *rng_within_site_p;
        orig_call call = orig_call(rng_within_site, is_gvcf);
        orig_calls_for_site[call] += n_calls;
    }
    return Status::OK();
}

// Called after joint genotyping of a unified site:
// n_no_calls gives the count of no-calls in the output joint call.
Status LossTracker::finalize_loss_for_site(int n_no_calls) noexcept {
    if (is_finalized)
        return Status::Invalid("calling finalize_loss_for_site when LossTracker is already finalized.");

    n_no_calls_total += n_no_calls;
    for (auto& kv : orig_calls_for_site) {
        int call_within_site_len = kv.first.pos.size();
        int n_orig_calls = kv.second;

        // Update total coverage of original calls implicated
        // in joint call for this site
        n_calls_total += n_orig_calls;
        n_bp_total += n_orig_calls * call_within_site_len;

        if (kv.first.is_gvcf) {
            n_gvcf_calls_total += n_orig_calls;
            n_gvcf_bp_total += n_orig_calls * call_within_site_len;
        }

        // Joint call has at least 1 missing call.
        if (n_no_calls) {
            // The expected behavior for computing n_calls_lost_for_site:
            //   If output joint call has 1 no-call:
            //      if n_orig_calls = 1 --> 0 lost calls
            //      if n_orig_calls = 2 --> 1 lost call
            //   If output joint call has 2 no-calls:
            //      if n_orig_calls = 1 --> 1 lost call
            //      if n_orig_calls = 2 --> 2 lost calls
            // It is not expected for n_orig_calls > 2 (this may happen if
            // multiple original records cover the same range after
            // intersecting with site). In this unlikely case, we compute
            // n_calls_lost_for_site as n_calls_lost multipled by
            // n_orig_calls divided by 2, rounded down to the nearest int

            int n_calls_lost_for_site = (n_orig_calls * n_no_calls) / 2;
            n_calls_lost += n_calls_lost_for_site;

            // call_within_site_len gives length of orig_call
            // restricted to the unified_site. Number of base pairs
            // (of original calls) lost is given by this length multipled
            // by n_calls_lost_for_site computed above
            n_bp_lost += call_within_site_len * n_calls_lost_for_site;

            // gvcf loss accounting
            if (kv.first.is_gvcf) {
                n_gvcf_calls_lost += n_calls_lost_for_site;
                n_gvcf_bp_lost += call_within_site_len * n_calls_lost_for_site;
            }
        }
    }

    // Clear map
    orig_calls_for_site.clear();
    is_finalized = true;

    return Status::OK();
}

// Returns the count variables packaged within loss_stats
Status LossTracker::get(loss_stats& ans) const noexcept {
    if (!is_finalized){
        return Status::Invalid("calling get on an unfinalized LossTracker");
    }

    ans.n_calls_total = n_calls_total;
    ans.n_bp_total = n_bp_total;
    ans.n_calls_lost = n_calls_lost;
    ans.n_no_calls_total = n_no_calls_total;
    ans.n_bp_lost = n_bp_lost;
    ans.n_gvcf_bp_lost = n_gvcf_bp_lost;
    ans.n_gvcf_calls_lost = n_gvcf_calls_lost;
    ans.n_gvcf_bp_total = n_gvcf_bp_total;
    ans.n_gvcf_calls_total = n_gvcf_calls_total;

    return Status::OK();
}


// Update the loss_stats data structure with call information for
// original calls associated with a unified site
static Status update_orig_calls_for_loss(const genotyper_config& cfg, const vector<shared_ptr<bcf1_t>>& records, int n_bcf_samples, const bcf_hdr_t* dataset_header, const map<int,int>& sample_mapping, LossTrackers& losses_for_site) {
    for (auto& record: records) {
        Status s;
        range rng(record);
        int *gt = nullptr, gtsz = 0;
        int nGT = bcf_get_genotypes(dataset_header, record.get(), &gt, &gtsz);
        if (!gt || nGT != 2*n_bcf_samples) return Status::Failure("genotyper::update_orig_calls_for_loss bcf_get_genotypes");

        for (int i = 0; i < n_bcf_samples; i++) {
            int sample_ind = sample_mapping.at(i);
            auto& loss = losses_for_site[sample_ind];

            int n_calls = !bcf_gt_is_missing(gt[i*2]) + !bcf_gt_is_missing(gt[i*2 + 1]);
            loss.add_call_for_site(rng, n_calls, is_gvcf_ref_record(cfg, record.get()));
        }

        free(gt);
    }
    return Status::OK();
}

// Update the loss_stats data sturcture with the joint call for
// the unified site and finalize the loss measures
static Status update_joint_call_loss(bcf1_t* record, int n_bcf_samples, const vector<one_call>& gt, LossTrackers& losses_for_site) {

    if(n_bcf_samples != losses_for_site.size()) {
        return Status::Failure("update_joint_call_loss: number of samples and bcf does not match");
    }
    range rng(record);
    Status s;

    for (int i = 0; i < n_bcf_samples; i++) {
        auto& loss = losses_for_site[i];

        int n_gt_missing = (bcf_gt_is_missing(gt[i*2].allele) + bcf_gt_is_missing(gt[i*2 + 1].allele));

        assert(n_gt_missing <= 2);
        // Lock down the loss associated with this unified_site
        loss.finalize_loss_for_site(n_gt_missing);
    }

    return Status::OK();
}

// Translate the hard-called genotypes from the bcf1_t into our genotype
// vector, based on a mapping from the bcf1_t sample indices into indices of
// the genotype vector. Calls on update_orig_calls_for_loss to register
// original calls in the bcf1_t for loss calculations
static Status translate_genotypes(const genotyper_config& cfg, const unified_site& site,
                                  const string& dataset, const bcf_hdr_t* dataset_header,
                                  bcf1_t* record, const map<int,int>& sample_mapping,
                                  vector<one_call>& genotypes, const vector<int>& min_ref_depth,
                                  LossTrackers& losses_for_site) {
    assert(genotypes.size() == 2*min_ref_depth.size());
    Status s;

    // map the BCF's alleles onto the unified alleles
    vector<int> allele_mapping(record->n_allele, -1);

    // reference allele maps if it contains the unified site
    range rng(record);
    if (rng.overlaps(site.pos)) {
        allele_mapping[0] = 0;
    }

    // map the bcf1_t alt alleles according to unification
    // checking for valid dna regex match
    for (int i = 1; i < record->n_allele; i++) {
        string al(record->d.allele[i]);
        if (regex_match(al, regex_dna)) {
            auto p = site.unification.find(allele(rng, al));
            if (p != site.unification.end()) {
                allele_mapping[i] = p->second;
            }
        }
    }

    // get the genotype calls
    int *gt = nullptr, gtsz = 0;
    int nGT = bcf_get_genotypes(dataset_header, record, &gt, &gtsz);
    int n_bcf_samples = bcf_hdr_nsamples(dataset_header);
    if (!gt || nGT != 2*n_bcf_samples) return Status::Failure("genotyper::translate_genotypes bcf_get_genotypes");
    assert(record->n_sample == bcf_hdr_nsamples(dataset_header));

    // Update the depth of coverage info
    unique_ptr<AlleleDepthHelper> depth;
    S(AlleleDepthHelper::Open(cfg, dataset, dataset_header, record, depth));

    // for each shared sample, record the genotype call.
    for (const auto& ij : sample_mapping) {
        assert(2*ij.first < nGT);
        assert(ij.second < min_ref_depth.size());

        #define fill_allele(ofs)                                           \
            if (gt[2*ij.first+ofs] != bcf_int32_vector_end &&              \
                !bcf_gt_is_missing(gt[2*ij.first+ofs])) {                  \
                auto al = bcf_gt_allele(gt[2*ij.first+ofs]);               \
                assert(al >= 0 && al < record->n_allele);                  \
                if (depth->sufficient(ij.first, al)                        \
                    && depth->sufficient_ref(min_ref_depth[ij.second])) {  \
                    if (allele_mapping[al] >= 0) {                         \
                        genotypes[2*ij.second+(ofs)] =                     \
                            one_call(bcf_gt_unphased(allele_mapping[al]),  \
                                     NoCallReason::N_A);                   \
                    } else {                                               \
                        genotypes[2*ij.second+(ofs)].RNC =                 \
                            NoCallReason::LostAllele;                      \
                    }                                                      \
                } else {                                                   \
                    genotypes[2*ij.second+(ofs)].RNC =                     \
                        NoCallReason::InsufficientDepth;                   \
                }                                                          \
            }
        fill_allele(0)
        fill_allele(1)
    }

    if (gt) {
        free(gt);
    }

    return Status::OK();
}

/// Given a set of records overlapping a site, we find the
/// "representative record(s)" for a given site. The representative
/// record(s) found are returned via the ans ptr. This function also
/// modifies 2 genotyper substrates: min_ref_depth and genotypes.
/// min_ref_depth tracks the minimum depth of ref-records for a given sample,
/// and will be set to -1 if there are no ref records processed for a sample
/// min_ref_depth is expected to be initialized to -1 for all samples.
/// Modifications to genotypes vector is explained below
/// The representative record is expected to be passed into translate_genotypes for genotyping
/// =====================================================
/// Pre conditions:
/// min_ref_depth should be initialized to -1 for all samples
/// ans should be empty when passed in
/// ======================================================
/// Expected behaviors:
/// No records given
///      ans returned as empty, no modification to genotypes and failed_ref_depth vectors
/// Records do not span the entire range of site
///      ans returned as empty, min_ref_depth updated accordingly, no change to genotypes
/// Records span entire range, and consist of all gvcf records
///      ans set as the first gvcf record passed, min_ref_depth updated accordingly, no change to genotypes
/// Records span entire range, and consist of exactly 1 vcf variant record
///      ans set as the 1 vcf variant record, min_ref_depth updated accordingly, no change to genotypes
/// Records span entire range, and consists of >1 vcf variant record
///      ans returned as empty, min_ref_depth updated accordingly, genotypes for all samples in this dataset are updated to be lost calls with RNC=LostAllele

Status find_rep_record(const genotyper_config& cfg, const unified_site& site,
                       const vector<shared_ptr<bcf1_t>>& records, const string& dataset,
                       const bcf_hdr_t* hdr, int bcf_nsamples,
                       const map<int, int>& sample_mapping, vector<one_call>& genotypes,
                       vector<int>& min_ref_depth, shared_ptr<bcf1_t>& ans) {


    if (ans) {
        return Status::Invalid("find_rep_record: ans is not empty when passed into the function.");
    }
    // Clear ans to be empty
    ans.reset();

    Status s;
    unique_ptr<AlleleDepthHelper> depth;

    // To track the total range covered by the records
    vector<range> record_rngs;
    for (auto& record: records) {
        range record_rng(record.get());
        assert(record_rng.overlaps(site.pos));
        record_rngs.push_back(record_rng);
        S(AlleleDepthHelper::Open(cfg, dataset, hdr, record.get(), depth));

        // Update min_ref_depth vector
        for (int i=0; i<bcf_nsamples; i++) {
            // TODO: Consider support for hom ref records
            int ref_depth = depth->get_gvcf_depth(i);
            int mapped_sample = sample_mapping.at(i);
            assert(mapped_sample < min_ref_depth.size());

            // relevant gvcf confidence record
             if (ref_depth >= 0) {
                if (min_ref_depth[mapped_sample] < 0) {
                    min_ref_depth[mapped_sample] = ref_depth;
                } else {
                    min_ref_depth[mapped_sample] = min(min_ref_depth[mapped_sample], ref_depth);
                }

            }
        }
    }

    // records do not span the site, return ans as empty
    if (!site.pos.spanned_by(record_rngs)) {
        // The RNCs default to MissingData, which is appropriate if records is
        // empty. If however there are some records (but not enough to cover
        // the site), update the RNCs to PartialData.
        if (!record_rngs.empty()) {
            for (int i = 0; i < bcf_nsamples; i++) {
                genotypes[sample_mapping.at(i)*2].RNC =
                    genotypes[sample_mapping.at(i)*2+1].RNC = NoCallReason::PartialData;
            }
        }
        return Status::OK();
    }

    // continue to find representative record if records span site
    vector<shared_ptr<bcf1_t>> non_gvcf_records;
    copy_if(records.begin(), records.end(), back_inserter(non_gvcf_records), [&cfg](auto& record) {return !is_gvcf_ref_record(cfg, record.get()); });

    if (non_gvcf_records.size() == 0) {
        // All records spanning site are gvcf, we can use any gvcf record
        // as the representative, and translate_genotype will interpret  site
        // as wholly REF
        ans = records[0];
    } else  if (non_gvcf_records.size() == 1) {
        // If there is only 1 non_gvcf (ie vcf variant) record, this record
        // should be used as a representative record for this site
        ans = non_gvcf_records[0];
    } else {
        // If there are multiple variant vcf records, we cannot effectively
        // handle genotyping due to phase assertions, so we return an empty
        // vector as the representative record.

        // Update genotype vector's RNC. Distinguish UnphasedVariants from
        // OverlappingVariants
        // TODO: replace quadratic loop (but the # records is probably too
        // small to matter)
        auto rnc = NoCallReason::UnphasedVariants;
        for (int i = 0; i < records.size() && rnc == NoCallReason::UnphasedVariants; i++) {
            range rng_i(records[i]);
            for (int j = i+1; j < records.size(); j++) {
                if (rng_i.overlaps(range(records[j]))) {
                    rnc = NoCallReason::OverlappingVariants;
                    break;
                }
            }
        }
        for (int i = 0; i < bcf_nsamples; i++) {
            genotypes[sample_mapping.at(i)*2].RNC = rnc;
            genotypes[sample_mapping.at(i)*2+1].RNC = rnc;
        }
    }

    return Status::OK();
}

Status genotype_site(const genotyper_config& cfg, MetadataCache& cache, BCFData& data, const unified_site& site,
                     const std::string& sampleset, const vector<string>& samples,
                     const bcf_hdr_t* hdr, shared_ptr<bcf1_t>& ans, consolidated_loss& losses_for_site,
                     atomic<bool>* ext_abort) {
	Status s;

    // Initialize a vector for the unified genotype calls for each sample,
    // starting with everything missing. We'll then loop through BCF records
    // overlapping this site and fill in the genotypes as we encounter them.
    vector<one_call> genotypes(2*samples.size());

    LossTrackers loss_trackers;
    for (const auto& sample : samples) {
        loss_trackers.push_back(LossTracker(site.pos));
    }

    shared_ptr<const set<string>> samples2, datasets;
    vector<unique_ptr<RangeBCFIterator>> iterators;
    S(data.sampleset_range(cache, sampleset, site.pos,
                           samples2, datasets, iterators));
    assert(samples.size() == samples2->size());

    // for each pertinent dataset
    for (const auto& dataset : *datasets) {
        if (ext_abort && *ext_abort) {
            return Status::Aborted();
        }

        // load BCF records overlapping the site by "merging" the iterators
        shared_ptr<const bcf_hdr_t> dataset_header;
        vector<vector<shared_ptr<bcf1_t>>> recordss(iterators.size());
        vector<shared_ptr<bcf1_t>> records;

        for (const auto& iter : iterators) {
            string this_dataset;
            vector<shared_ptr<bcf1_t>> these_records;
            S(iter->next(this_dataset, dataset_header, these_records));
            if (dataset != this_dataset) {
                return Status::Failure("genotype_site: iterator returned unexpected dataset",
                                       this_dataset + " instead of " + dataset);
            }
            records.insert(records.end(), these_records.begin(), these_records.end());
        }

        assert(is_sorted(records.begin(), records.end(),
                         [] (shared_ptr<bcf1_t>& p1, shared_ptr<bcf1_t>& p2) {
                            return range(p1) < range(p2);
                         }));

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

        update_orig_calls_for_loss(cfg, records, bcf_nsamples, dataset_header.get(), sample_mapping, loss_trackers);

        vector<int> min_ref_depth(samples.size(), -1);
        shared_ptr<bcf1_t> rep_record;
        S(find_rep_record(cfg, site, records, dataset, dataset_header.get(), bcf_nsamples, sample_mapping, genotypes, min_ref_depth, rep_record));

        // non-empty rep_record
        if (rep_record) {
            S(translate_genotypes(cfg, site, dataset, dataset_header.get(), rep_record.get(),
                                  sample_mapping, genotypes, min_ref_depth, loss_trackers));
        }
    }

    // Clean up emission order of alleles
    for(size_t i=0; i < samples.size(); i++) {
        if(genotypes[2*i] > genotypes[2*i + 1]) {
            swap(genotypes[2*i], genotypes[2*i+1]);
        }
    }
    // Create the destination BCF record for this site.
    ans = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    ans->rid = site.pos.rid;
    ans->pos = site.pos.beg;
    ans->rlen = site.pos.end - site.pos.beg;
    ans->qual = 0;

    // alleles
    vector<const char*> c_alleles;
    for (const auto& allele : site.alleles) {
        c_alleles.push_back(allele.c_str());
    }
    if (bcf_update_alleles(hdr, ans.get(), c_alleles.data(), c_alleles.size()) != 0) {
        return Status::Failure("bcf_update_alleles");
    }

    // GT
    vector<int32_t> gt;
    for (const auto& c : genotypes) {
        gt.push_back(c.allele);
    }
    assert(gt.size() == genotypes.size());
    if (bcf_update_genotypes(hdr, ans.get(), gt.data(), gt.size()) != 0) {
        return Status::Failure("bcf_update_genotypes");
    }

    // RNC
    vector<const char*> rnc;
    for (const auto& c : genotypes) {
        char* v = "M";
        switch (c.RNC) {
            case NoCallReason::N_A:
                v = ".";
                break;
            case NoCallReason::PartialData:
                v = "P";
                break;
            case NoCallReason::LostAllele:
                v = "L";
                break;
            case NoCallReason::InsufficientDepth:
                v = "D";
                break;
            case NoCallReason::UnphasedVariants:
                v = "U";
                break;
            case NoCallReason::OverlappingVariants:
                v = "O";
                break;
            default:
                assert(c.RNC == NoCallReason::MissingData);
        }
        rnc.push_back(v);
    }
    assert (gt.size() == rnc.size());
    if (bcf_update_format_string(hdr, ans.get(), "RNC", rnc.data(), rnc.size()) != 0) {
        return Status::Failure("bcf_update_format_string RNC");
    }

    // Finalize loss statistics for this site
    S(update_joint_call_loss(ans.get(), bcf_hdr_nsamples(hdr), genotypes, loss_trackers));
    // Package consolidated_loss for this site and merge into losses_for_site
    // to be returned to parent caller
    consolidated_loss losses;
    for (int i = 0; i < loss_trackers.size(); i++) {
        auto& tracker = loss_trackers[i];
        auto& sample_name = samples[i];
        loss_stats loss;
        S(tracker.get(loss));
        losses.insert(make_pair(sample_name,loss));
    }
    merge_loss_stats(losses, losses_for_site);
    return Status::OK();
}

}
