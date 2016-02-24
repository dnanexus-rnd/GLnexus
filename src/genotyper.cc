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
    return record->n_allele == 2 && strcmp(record->d.allele[1], cfg.ref_symbolic_allele.c_str()) == 0;
}

/// Detect an idiosyncratic class of records from certain HaplotypeCaller
/// versions which have QUAL == 0.0 and 0/0 genotype calls...we treat these as
/// "pseudo" reference confidence records.
static bool is_pseudo_ref_record(const bcf_hdr_t* hdr, bcf1_t* record) {
    if (record->qual != 0.0) {
        return false;
    }
    int *gt = nullptr, gtsz = 0;
    int nGT = bcf_get_genotypes(hdr, record, &gt, &gtsz);
    for (int i = 0; i < record->n_sample; i++) {
        assert(i*2+1<nGT);
        if (bcf_gt_is_missing(gt[i*2]) || bcf_gt_allele(gt[i*2]) != 0 ||
            bcf_gt_is_missing(gt[i*2+1]) || bcf_gt_allele(gt[i*2+1]) != 0) {
            free(gt);
            return false;
        }
    }
    free(gt);
    return true;
}

class IFormatFieldHelper {

public:
    // basic information for the retained field
    retained_format_field field_info;


    // Expected number of samples in output bcf record
    int n_samples;

    // Expected number of values per sample in the **output** (ie unified site)
    // Note this may differ from the count of input for RetainedFieldNumber::ALLELES
    // and RetainedFieldNumber::ALT since the number of alleles may differ in input and
    // output
    int count;

     IFormatFieldHelper(const retained_format_field field_info_, int n_samples_, int count_) : field_info(field_info_), n_samples(n_samples_), count(count_) {}

     IFormatFieldHelper() = default;

     virtual Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header, bcf1_t* record,
                                    const map<int, int>& sample_mapping, const vector<int> allele_mapping) = 0;

     virtual Status update_record_format(const bcf_hdr_t* hdr, bcf1_t* record) = 0;

     virtual ~IFormatFieldHelper() = default;

};

template <class T>
class FormatFieldHelper : public IFormatFieldHelper {

    // A vector of length n_samples, where each element
    // is a vector consisting of format field value from
    // record(s) related to the given sample
    vector<vector<T>> format_v;

    // Combination function to handle combining multiple format values
    // from multiple records
    T (*combine_f) (vector<T>);

    static T max_element_wrapper(vector<T> v) {
        return (*max_element(v.begin(), v.end()));
    }

    static T min_element_wrapper(vector<T> v) {
        return (*min_element(v.begin(), v.end()));
    }

    // Overloaded wrapper function to call bcf_get_format of the correct
    // format field type
    static int bcf_get_format_wrapper(const bcf_hdr_t* dataset_header, bcf1_t* record, const char* field_name, int32_t** v, int* vsz) {
        return bcf_get_format_int32(dataset_header, record, field_name, v, vsz);
    }
    static int bcf_get_format_wrapper(const bcf_hdr_t* dataset_header, bcf1_t* record, const char* field_name, float** v, int* vsz) {
        return bcf_get_format_float(dataset_header, record, field_name, v, vsz);
    }

    Status combine_format_data(vector<T>& ans) {
        if (field_info.default_to_zero) {
            for (auto& v : format_v) {
                if (v.empty()) {
                    v.push_back(0);
                }
            }
        }

        if( any_of(format_v.begin(), format_v.end(), [](vector<T> v){return v.empty();})) {
            return Status::Invalid("genotyper: one or more sample has missing FORMAT field for the intended output field", field_info.name);
        }
        assert(format_v.size() == n_samples * count);

        for (auto& format_one : format_v) {
            ans.push_back(combine_f(format_one));
        }

        return Status::OK();
    }

    /// For a given record, give the number of expected values
    /// per sample, based on the format type
    int expected_n_val_per_sample(const bcf1_t* record) {

        // RetainedFieldNumber::BASIC
        int expected_count = count;
        if (field_info.number == RetainedFieldNumber::ALT) {
            expected_count = record->n_allele - 1;
        } else if (field_info.number == RetainedFieldNumber::ALLELES) {
            expected_count = record->n_allele;
        } else if (field_info.number == RetainedFieldNumber::GENOTYPE) {
            expected_count = (record->n_allele + 1) * record->n_allele / 2;
        }
        return expected_count;
    }

    /// Given a format field for an input sample (unmapped_i), and
    /// the unmapped_j-th value for this sample, find the corresponding
    /// index of the this format field in the output.
    /// Returns a negative value if this value cannot be mapped to the output
    /// (e.g. allele-specific info for a trimmed allele), and raises error
    /// if the sample cannot be mapped
    int get_out_ind_of_value(int unmapped_i, int unmapped_j,
                             const map<int, int>& sample_mapping,
                             const vector<int> allele_mapping) {
        int mapped_i = sample_mapping.at(unmapped_i);

        // Sample should always be mappable
        assert (mapped_i >= 0);

        switch (field_info.number){
            case RetainedFieldNumber::ALT:
            case RetainedFieldNumber::ALLELES:
            {
                // Fall through for both cases that require allele_mapping
                int mapped_j = allele_mapping[unmapped_j];
                if (mapped_j < 0) {
                    // Allele is not mappable (trimmed or is a gvcf record)
                    return -1;
                } else {
                    return mapped_i * count + mapped_j;
                }
                break;
            }
            case RetainedFieldNumber::GENOTYPE:
            {
                // TODO: Generate mapping for PL field (Number=G)
                return -1;
            }
            default:
            {
                // RetainedFieldNumber::BASIC case
                return mapped_i * count + unmapped_j;
            }
        }
    }

    // This code is a duplication of add_record_data, for the most part,
    // I'm leaving it as duplicated first, for future considerations on making this more compact
    // after we have more of a handle on how to deal with the AD field
    Status handle_AD_field(const string& dataset, const bcf_hdr_t* dataset_header, bcf1_t* record,
                           const map<int, int>& sample_mapping, const vector<int> allele_mapping) {

        // Fill AD using either DP or MIN_DP, preferentially use DP, if available
        vector<string> depth_fields = {"DP", "MIN_DP"};
        bool found = false;

        cout << "Handling AD case for " << record->rid << " " << record->pos << endl;
        // Expect only 1 value (corresponding the REF record)
        // This conveniently help us with the allele-mapping step, since the
        // first allele is guaranteed to be the REF allele
        int n_val_per_sample = 1;

        for (auto& field_name : depth_fields) {
            if (found) break;

            T *v = nullptr;
            int vsz = 0;

            // rv is the number of values written
            int rv = FormatFieldHelper::bcf_get_format_wrapper(dataset_header, record, field_name.c_str(), &v, &vsz);

            // raise error if there's a failed get due to type mismatch
            if (rv == -2) {
                if (v) free(v);
                ostringstream errmsg;
                errmsg << dataset << " " << range(record).str() << " (" << field_name << ")";
                return Status::Invalid("genotyper: getting format field errored with type mismatch", errmsg.str());
            }
            // don't raise error if get failed due to tag missing in record or
            // in vcf header; continue to look with other possible field names

            if (rv >= 0) {
                found = true;
                if (rv != record->n_sample * n_val_per_sample) {
                // For this field, we expect n_val_per_sample values per sample
                    if (v) free(v);
                    ostringstream errmsg;
                    errmsg << dataset << " " << range(record).str() << "(" << field_name << ")";
                    return Status::Invalid("genotyper: unexpected result when fetching record FORMAT field", errmsg.str());
                } // close rv != record->n_sample * count

                for (int i=0; i<record->n_sample; i++) {
                    for (int j=0; j<n_val_per_sample; j++) {

                        int in_ind = i * n_val_per_sample + j;
                        int out_ind = get_out_ind_of_value(i, j, sample_mapping, allele_mapping);

                        if (out_ind < 0) {
                           continue;
                        }

                        assert(out_ind < format_v.size());
                        assert(in_ind < vsz);
                        format_v[out_ind].push_back(v[in_ind]);
                    } // close for j loop
                } // close for i loop
            } // close rv >= 0
            free(v);
        }

        return Status::OK();
    }
public:

    FormatFieldHelper(const retained_format_field field_info_, int n_samples_, int count_) : IFormatFieldHelper(field_info_, n_samples_, count_) {

        switch (field_info.combi_method) {
            case FieldCombinationMethod::MIN:
                combine_f = FormatFieldHelper::min_element_wrapper;
                break;
            case FieldCombinationMethod::MAX:
                combine_f = FormatFieldHelper::max_element_wrapper;
                break;
        }

        for (int i=0; i<n_samples_ * count_; i++){
            // We keep 1 vector for each eventual output value
            // So if every sample contains X values (in vcf specification, Number=X), then we keep X vectors per sample.
            format_v.push_back(vector<T>());
        }
    }

    virtual ~FormatFieldHelper() = default;


    Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header, bcf1_t* record,
                           const map<int, int>& sample_mapping, const vector<int> allele_mapping) {

        bool found = false;
        int n_val_per_sample = expected_n_val_per_sample(record);
        for (auto& field_name : field_info.orig_names) {
            if (found) break;
            T *v = nullptr;
            int vsz = 0;

            // rv is the number of values written
            int rv = FormatFieldHelper::bcf_get_format_wrapper(dataset_header, record, field_name.c_str(), &v, &vsz);

            // raise error if there's a failed get due to type mismatch
            if (rv == -2) {
                if (v) free(v);
                ostringstream errmsg;
                errmsg << dataset << " " << range(record).str() << " (" << field_name << ")";
                return Status::Invalid("genotyper: getting format field errored with type mismatch", errmsg.str());
            }
            // don't raise error if get failed due to tag missing in record or
            // in vcf header; continue to look with other possible field names

            if (rv >= 0) {
                found = true;
                if (rv != record->n_sample * n_val_per_sample) {
                // For this field, we expect n_val_per_sample values per sample
                    if (v) free(v);
                    ostringstream errmsg;
                    errmsg << dataset << " " << range(record).str() << "(" << field_name << ")";
                    return Status::Invalid("genotyper: unexpected result when fetching record FORMAT field", errmsg.str());
                } // close rv != record->n_sample * count

                for (int i=0; i<record->n_sample; i++) {
                    for (int j=0; j<n_val_per_sample; j++) {

                        int in_ind = i * n_val_per_sample + j;
                        int out_ind = get_out_ind_of_value(i, j, sample_mapping, allele_mapping);

                        if (out_ind < 0) {
                           continue;
                        }

                        assert(out_ind < format_v.size());
                        assert(in_ind < vsz);
                        format_v[out_ind].push_back(v[in_ind]);
                    } // close for j loop
                } // close for i loop
            } // close rv >= 0
            free(v);
        }

        if (!found) {

            // Special case handling for AD field to convert ref
            // DP/MIN_DP to repopulate AD field
            if (field_info.name == "AD") {
                return handle_AD_field(dataset, dataset_header, record, sample_mapping, allele_mapping);
            }
            // TODO: The specified field could not be found from the input record, this is "OK"
            // because we might be able to find some other record that has the requisite value
            // Populating by default value is handled in the update_record_format step
            return Status::OK();
        }
        return Status::OK();
    }

    Status update_record_format(const bcf_hdr_t* hdr, bcf1_t* record) {
        Status s;
        vector<T> ans;
        s = combine_format_data(ans);

        if (!s.ok()) {
            // TODO: Combine function failed; currently ignore error, do not
            // update this field, and move on. More intelligent error handling?
            return Status::OK();
        }
        int retval  = 0;
        switch (field_info.type) {
            case RetainedFieldType::INT:
                retval = bcf_update_format_int32(hdr, record, field_info.name.c_str(), ans.data(), n_samples * count);
                break;
            case RetainedFieldType::FLOAT:
                retval = bcf_update_format_float(hdr, record, field_info.name.c_str(), ans.data(), n_samples * count);
                break;
            default:
                return Status::Invalid("genotyper: Unexpected RetainedFieldType when executing update_record_format.");
        }
        if (retval != 0) {
            return Status::Failure("genotyper: failed to update record format when executing update_record_format.");
        }
        return Status::OK();
    }


};


// Helper class for keeping track of the per-allele depth of coverage info in
// a bcf1_t record. There are a couple different cases to handle, depending on
// whether we're looking at a gVCF reference confidence record or a "regular"
// VCF record.
class AlleleDepthHelper {
    Status s_;
    const genotyper_config& cfg_;
    size_t n_sample_ = 0, n_allele_ = 0;
    bool is_g_ = false;
    int32_t *v_ = nullptr;
    int vsz_ = 0;

public:

    // The AlleleDepthHelper is constructed into an undefined state. Load()
    // must be invoked, successfully, before it can be used.
    AlleleDepthHelper(const genotyper_config& cfg)
        : cfg_(cfg)
        {}

    ~AlleleDepthHelper() {
        if (v_) free(v_);
    }

    // The helper can be reused for multiple records by calling Load()
    // repeatedly. This will be slightly more efficient than using a new
    // helper for each record.
    Status Load(const string& dataset, const bcf_hdr_t* dataset_header, bcf1_t* record) {
        n_sample_ = record->n_sample;
        n_allele_ = record->n_allele;
        // is this a gVCF reference confidence record?
        is_g_ = is_gvcf_ref_record(cfg_, record);

        if (is_g_) {
            // if so, look for the MIN_DP FORMAT field (or equivalent)
            int nv = bcf_get_format_int32(dataset_header, record, cfg_.ref_dp_format.c_str(),
                                          &v_, &vsz_);

            if (nv != record->n_sample) {
                ostringstream errmsg;
                errmsg << dataset << " " << range(record).str() << " (" << cfg_.ref_dp_format << ")";
                return Status::Invalid("genotyper: gVCF reference depth FORMAT field is missing or malformed", errmsg.str());
            }
        } else {
            // this is a regular VCF record, so look for the AD FORMAT field (or equivalent)
            int nv = bcf_get_format_int32(dataset_header, record, cfg_.allele_dp_format.c_str(),
                                          &v_, &vsz_);

            if (nv == -1) {
                // We allow the AD field to not exist mainly as a (poor)
                // workaround for some of our test case gVCFs not having it...
                size_t sz = record->n_sample * record->n_allele * sizeof(int32_t);
                if (vsz_ < sz) {
                    v_ = (int32_t*) realloc(v_, sz);
                    vsz_ = sz;
                }
                memset(v_, 0, sz);
            } else if (nv != record->n_sample * record->n_allele) {
                ostringstream errmsg;
                errmsg << dataset << " " << range(record).str() << " (" << cfg_.allele_dp_format << ")";
                return Status::Invalid("genotyper: VCF allele depth FORMAT field is malformed", errmsg.str());
            }
        }
        return Status::OK();
    }

    // The behavior of remaining methods is undefined until Load() has been
    // invoked successfully

    // Is there sufficient coverage for sample i, allele j?
    bool sufficient(unsigned sample, unsigned allele) {
        if (!cfg_.required_dp) {
            return true;
        }
        assert(sample < n_sample_);
        assert(allele < n_allele_);
        if (is_g_) {
            // the MIN_DP array has just one integer per sample
            if (allele == 0) {
                return v_[sample] >= cfg_.required_dp;
            } else {
                return false;
            }
        } else {
            // the AD array has one integer per allele per sample
            return v_[sample*n_allele_+allele] >= cfg_.required_dp;
        }
    }

    bool is_gvcf_ref() { return is_g_; }

    // Returns the depth for the reference allele for sample i
    int get_ref_depth(unsigned sample) {
        assert(sample < n_sample_);
        if (is_g_) {
            return v_[sample];
        } else {
            return v_[sample*n_allele_+0];
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

// Helper: given REF and ALT DNA, determine if the ALT represents a deletion
// with respect to REF. Left-alignment is assumed and reference padding on the
// left is tolerated.
inline bool is_deletion(const string& ref, const string& alt) {
    if (alt.size() >= ref.size()) {
        return false;
    }
    for (int j = 0; j < alt.size(); j++) {
        if (alt[j] != ref[j]) {
            // some kind of complex edit where the ALT allele is both shorter
            // and with differing basis in the prefix...
            return false;
        }
    }
    return true;
}

/// A helper function to update min_ref_depth based on several reference
/// confidence records. min_ref_depth[j] is the minimum depth of reference
/// coverage seen for sample j across the reference confidence records, and
/// should be initialized to -1 before any reference confidence records are
/// seen.
static Status update_min_ref_depth(const string& dataset, const bcf_hdr_t* dataset_header,
                                   int bcf_nsamples, const map<int,int>& sample_mapping,
                                   const vector<shared_ptr<bcf1_t>>& ref_records,
                                   AlleleDepthHelper& depth,
                                   vector<int>& min_ref_depth) {
    Status s;
    for (auto& ref_record : ref_records) {
        S(depth.Load(dataset, dataset_header, ref_record.get()));

        for (int i=0; i<bcf_nsamples; i++) {
            int mapped_sample = sample_mapping.at(i);
            assert(mapped_sample < min_ref_depth.size());

            if (min_ref_depth[mapped_sample] < 0) {
                min_ref_depth[mapped_sample] = depth.get_ref_depth(i);
            } else {
                min_ref_depth[mapped_sample] = min(min_ref_depth[mapped_sample],
                                                   depth.get_ref_depth(i));
            }
        }
    }
    return Status::OK();
}

/// Given a unified site and the set of gVCF records overlapping it in some
/// dataset, separate the variant records from surrounding reference
/// confidence records. Ideally and often there's either zero or one variant
/// records -- zero if the dataset exhibits no variation at the site, and one
/// if it does. Unfortunately, there are a number of circumstances under which
/// some variant callers produce multiple overlapping records. The variant
/// records should all share at least one reference position in common.
///
/// The variant records, if any, are returned through variant_records. It is
/// then the job of translate_genotypes to figure out what to do with the
/// cluster of variant records.
///
/// This function also modifies min_ref_depth which tracks the minimum depth
/// of ref-records for a given sample. min_ref_depth is expected to be
/// initialized to -1 for all samples, and the entry for a sample will be
/// modified only if one or more ref records are present.
///
/// =====================================================
/// Pre conditions:
/// min_ref_depth should be initialized to -1 for all samples
/// ======================================================
/// Post conditions:
/// No records given
///      variant_records empty, min_ref_depth unchanged, rnc = MissingData
/// Records do not span the entire range of site
///      variant_records empty, min_ref_depth updated accordingly, rnc = PartialData
/// Records span entire range, and consist of all reference confidence records
///      variant_records empty, min_ref_depth updated accordingly, rnc = N_A
/// Records span entire range, and include one or more variant records which
/// all share at least one reference position
///      variant_records filled in, min_ref_depth updated accordingly, rnc = N_A
/// Records span entire range, and include multiple variant records which
/// don't all share at least one reference position
///      variant_records empty, min_ref_depth updated accordingly, rnc = UnphasedVariants
Status find_variant_records(const genotyper_config& cfg, const unified_site& site,
                            const string& dataset, const bcf_hdr_t* hdr, int bcf_nsamples,
                            const map<int, int>& sample_mapping,
                            const vector<shared_ptr<bcf1_t>>& records,
                            AlleleDepthHelper& depth,
                            NoCallReason& rnc,
                            vector<int>& min_ref_depth,
                            vector<shared_ptr<bcf1_t>>& variant_records) {
    // initialize outputs
    rnc = NoCallReason::MissingData;
    variant_records.clear();

    Status s;

    // collect the ranges covered by the records
    vector<range> record_rngs;
    record_rngs.reserve(records.size());
    for (auto& record: records) {
        range record_rng(record.get());
        assert(record_rng.overlaps(site.pos));
        record_rngs.push_back(record_rng);
        assert(bcf_nsamples == record->n_sample);
    }

    // check the records span the site, otherwise we need to produce
    // PartialData no-calls
    if (!site.pos.spanned_by(record_rngs)) {
        if (!record_rngs.empty()) {
            rnc = NoCallReason::PartialData;
        }
        return Status::OK();
    }

    // now that we know the records span the site, partition the variant
    // record(s) and the surrounding reference confidence records.
    vector<shared_ptr<bcf1_t>> ref_records;
    ref_records.reserve(records.size());
    for (auto& record : records) {
        (is_gvcf_ref_record(cfg, record.get()) ? ref_records : variant_records).push_back(record);
    }

    // compute min_ref_depth across the reference confidence records
    S(update_min_ref_depth(dataset, hdr, bcf_nsamples, sample_mapping,
                           ref_records, depth, min_ref_depth));

    // if there are multiple variant records, check if they all share at least
    // one reference position.
    if (variant_records.size() > 1) {
        range intersection(variant_records[0]);
        for (auto& record : variant_records) {
            range record_rng(record);
            assert(record_rng.rid == intersection.rid);
            intersection.beg = max(record_rng.beg, intersection.beg);
            intersection.end = min(record_rng.end, intersection.end);
        }
        if (intersection.beg >= intersection.end) {
            // if not, then we've got to bug out with UnphasedVariants
            variant_records.clear();
            rnc = NoCallReason::UnphasedVariants;
            return Status::OK();
        }
    }

    // Success...
    rnc = NoCallReason::N_A;
    return Status::OK();
}

static Status find_allele_mapping(const unified_site& site, const bcf1_t *record,
                                  vector<int>& allele_mapping, vector<bool>& deletion_allele) {
    range rng(record);
    assert(rng.overlaps(site.pos));
    allele_mapping[0] = 0;

    // map the bcf1_t alt alleles according to unification
    // checking for valid dna regex match
    string ref_al(record->d.allele[0]);
    for (int i = 1; i < record->n_allele; i++) {
        string al(record->d.allele[i]);
        if (regex_match(al, regex_dna)) {
            auto p = site.unification.find(allele(rng, al));
            if (p != site.unification.end()) {
                allele_mapping[i] = p->second;
            }
        }
        if (al.size() < rng.size() && rng.size() == ref_al.size()) {
            deletion_allele[i] = is_deletion(ref_al, al);
        }
    }

    return Status::OK();
}

/// Based on the cluster of variant records and min_ref_depth produced by
/// find_variant_records, fill genotypes for this dataset's samples with
/// appropriate calls (currently by translation of the input hard-calls).
/// Updates genotypes and losses_for_site, and may modify min_ref_depth.
static Status translate_genotypes(const genotyper_config& cfg, const unified_site& site,
                                  const string& dataset, const bcf_hdr_t* dataset_header,
                                  int bcf_nsamples, const map<int,int>& sample_mapping,
                                  const vector<shared_ptr<bcf1_t>>& records,
                                  AlleleDepthHelper& depth,
                                  vector<int>& min_ref_depth,
                                  vector<one_call>& genotypes,
                                  LossTrackers& losses_for_site) {
    assert(genotypes.size() == 2*min_ref_depth.size());
    Status s;

    // right now we can only actually deal with one variant record...
    bcf1_t *record = nullptr;

    if (records.size() == 1) {
        record = records[0].get();
    } else if (records.size() > 1) {
        // if we're given multiple overlapping variant records, attempt to
        // reclassify all, or all but one, of them as "pseudo" ref records
        vector<shared_ptr<bcf1_t>> pseudo_ref_records;
        for (auto& a_record : records) {
            assert(!is_gvcf_ref_record(cfg, a_record.get()));
            if (is_pseudo_ref_record(dataset_header, a_record.get())) {
                pseudo_ref_records.push_back(a_record);
            } else if (!record) {
                record = a_record.get();
            } else {
                // we've got an OverlappingVariants problem...
                for (int i = 0; i < bcf_nsamples; i++) {
                    genotypes[sample_mapping.at(i)*2].RNC =
                        genotypes[sample_mapping.at(i)*2+1].RNC =
                            NoCallReason::OverlappingVariants;
                }

                return Status::OK();
            }
        }

        assert(pseudo_ref_records.size());
        // update min_ref_depth with these pseudo ref records
        S(update_min_ref_depth(dataset, dataset_header, bcf_nsamples, sample_mapping,
                               pseudo_ref_records, depth, min_ref_depth));
    }

    if (!record) {
        // no variation represented in this dataset; make homozygous ref calls
        // if min_ref_depth indicates sufficient coverage for this sample
        for (const auto& ij : sample_mapping) {
            assert(ij.second < min_ref_depth.size());
            int rd = min_ref_depth[ij.second];

            if (rd >= cfg.required_dp) {
                genotypes[2*ij.second] =
                    genotypes[2*ij.second+1] =
                        one_call(bcf_gt_unphased(0), NoCallReason::N_A);
            } else {
                genotypes[2*ij.second].RNC =
                    genotypes[2*ij.second+1].RNC = NoCallReason::InsufficientDepth;
            }
        }
        return Status::OK();
    }

    // Now, translating genotypes from one variant BCF record.

    // map the BCF's alleles onto the unified alleles
    vector<int> allele_mapping(record->n_allele, -1);
    vector<bool> deletion_allele(record->n_allele, false); // which alleles are deletions

    S(find_allele_mapping(site, record, allele_mapping, deletion_allele))

    // range rng(record);
    // assert(rng.overlaps(site.pos));
    // allele_mapping[0] = 0;

    // // map the bcf1_t alt alleles according to unification
    // // checking for valid dna regex match
    // string ref_al(record->d.allele[0]);
    // for (int i = 1; i < record->n_allele; i++) {
    //     string al(record->d.allele[i]);
    //     if (regex_match(al, regex_dna)) {
    //         auto p = site.unification.find(allele(rng, al));
    //         if (p != site.unification.end()) {
    //             allele_mapping[i] = p->second;
    //         }
    //     }
    //     if (al.size() < rng.size() && rng.size() == ref_al.size()) {
    //         deletion_allele[i] = is_deletion(ref_al, al);
    //     }
    // }

    // get the genotype calls
    int *gt = nullptr, gtsz = 0;
    int nGT = bcf_get_genotypes(dataset_header, record, &gt, &gtsz);
    int n_bcf_samples = bcf_hdr_nsamples(dataset_header);
    if (!gt || nGT != 2*n_bcf_samples) return Status::Failure("genotyper::translate_genotypes bcf_get_genotypes");
    assert(record->n_sample == bcf_hdr_nsamples(dataset_header));

    S(depth.Load(dataset, dataset_header, record));

    // for each shared sample, record the genotype call.
    for (const auto& ij : sample_mapping) {
        assert(2*ij.first < nGT);
        assert(ij.second < min_ref_depth.size());

        #define fill_allele(ofs)                                           \
            if (gt[2*ij.first+ofs] != bcf_int32_vector_end &&              \
                !bcf_gt_is_missing(gt[2*ij.first+(ofs)])) {                \
                auto al = bcf_gt_allele(gt[2*ij.first+(ofs)]);             \
                assert(al >= 0 && al < record->n_allele);                  \
                int rd = min_ref_depth[ij.second];                         \
                if (depth.sufficient(ij.first, al)                         \
                    && (rd < 0 || rd >= cfg.required_dp)) {                \
                    if (allele_mapping[al] >= 0) {                         \
                        genotypes[2*ij.second+(ofs)] =                     \
                            one_call(bcf_gt_unphased(allele_mapping[al]),  \
                                     NoCallReason::N_A);                   \
                    } else {                                               \
                        genotypes[2*ij.second+(ofs)].RNC =                 \
                            deletion_allele[al]                            \
                                ? NoCallReason::LostDeletion               \
                                : NoCallReason::LostAllele;                \
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

    // Setup format field helpers
    vector<unique_ptr<IFormatFieldHelper>> format_helpers;
    for (const auto& format_field_info : cfg.liftover_fields) {
        int format_sz = samples.size();
        int count = -1;
        if (format_field_info.number == RetainedFieldNumber::BASIC) {
            count = format_field_info.count;
        } else if (format_field_info.number == RetainedFieldNumber::ALT) {
            // site.alleles.size() gives # alleles incl. REF
            count = (site.alleles.size() - 1);
        } else if (format_field_info.number == RetainedFieldNumber::ALLELES) {
            count = (site.alleles.size());
        }
        switch (format_field_info.type) {
            case RetainedFieldType::INT:
            {
                format_helpers.push_back(unique_ptr<IFormatFieldHelper>(new FormatFieldHelper<int32_t>(format_field_info, samples.size(), count)));
                break;
            }
            case RetainedFieldType::FLOAT:
            {
                format_helpers.push_back(unique_ptr<IFormatFieldHelper>(new FormatFieldHelper<float>(format_field_info, samples.size(), count)));
                break;
            }
        }
    }

    shared_ptr<const set<string>> samples2, datasets;
    vector<unique_ptr<RangeBCFIterator>> iterators;
    S(data.sampleset_range(cache, sampleset, site.pos,
                           samples2, datasets, iterators));
    assert(samples.size() == samples2->size());

    AlleleDepthHelper adh(cfg);

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

        // update loss trackers
        update_orig_calls_for_loss(cfg, records, bcf_nsamples, dataset_header.get(),
                                   sample_mapping, loss_trackers);

        // Update format fields
        for (auto& record : records) {
            vector<int> allele_mapping(record->n_allele, -1);
            vector<bool> deletion_allele(record->n_allele, false); // which alleles are deletions

            S(find_allele_mapping(site, record.get(), allele_mapping, deletion_allele))

            for (auto& format_helper : format_helpers) {
                format_helper->add_record_data(dataset, dataset_header.get(), record.get(),
                                               sample_mapping, allele_mapping);
            }
        }

        // find the variant records and process the surrounding reference
        // confidence records
        vector<int> min_ref_depth(samples.size(), -1);
        vector<shared_ptr<bcf1_t>> variant_records;
        NoCallReason rnc = NoCallReason::MissingData;
        S(find_variant_records(cfg, site, dataset, dataset_header.get(), bcf_nsamples,
                               sample_mapping, records, adh, rnc, min_ref_depth, variant_records));

        if (rnc != NoCallReason::N_A) {
            // no call for the samples in this dataset (several possible
            // reasons)
            for (int i = 0; i < bcf_nsamples; i++) {
                genotypes[sample_mapping.at(i)*2].RNC =
                    genotypes[sample_mapping.at(i)*2+1].RNC = rnc;
            }
        } else {
            // make genotype calls for the samples in this dataset
            S(translate_genotypes(cfg, site, dataset, dataset_header.get(), bcf_nsamples,
                                  sample_mapping, variant_records, adh, min_ref_depth,
                                  genotypes, loss_trackers));
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

    // Lifted-over FORMAT fields (non-genotype based)
    // TODO: It seems like the update command screw up the
    // alleles (REF, ALT) fields, quick fix: update the format
    // before updating the alleles, but KIV bug fix.
    for (auto& format_helper : format_helpers) {
        S(format_helper->update_record_format(hdr, ans.get()));
    }

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
        char* v = (char*) "M";
        #define RNC_CASE(reason,code) case NoCallReason::reason: v = (char*) code ; break;
        switch (c.RNC) {
            RNC_CASE(N_A,".")
            RNC_CASE(PartialData,"P")
            RNC_CASE(InsufficientDepth,"D")
            RNC_CASE(LostDeletion,"-")
            RNC_CASE(LostAllele,"L")
            RNC_CASE(UnphasedVariants,"U")
            RNC_CASE(OverlappingVariants,"O")
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
