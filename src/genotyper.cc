#include <assert.h>
#include <algorithm>
#include "genotyper.h"
#include "diploid.h"

using namespace std;

namespace GLnexus {

// Helper: determine if all genotypes in the record are 0/0
static inline bool is_homozygous_ref(const bcf_hdr_t* hdr, bcf1_t* record) {
    htsvecbox<int> gt;
    int nGT = bcf_get_genotypes(hdr, record, &gt.v, &gt.capacity);
    for (int i = 0; i < record->n_sample; i++) {
        assert(i*2+1<nGT);
        if (bcf_gt_is_missing(gt[i*2]) || bcf_gt_allele(gt[i*2]) != 0 ||
            bcf_gt_is_missing(gt[i*2+1]) || bcf_gt_allele(gt[i*2+1]) != 0) {
            return false;
        }
    }
    return true;
}

// Helper: given REF and ALT DNA, determine if the ALT represents a deletion
// with respect to REF. Left-alignment is assumed and reference padding on the
// left is tolerated.
static inline bool is_deletion(const string& ref, const string& alt) {
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

// Pre-process a bcf1_t record to cache some useful info that we'll use repeatedly
 Status preprocess_record(const unified_site& site, const bcf_hdr_t* hdr, const shared_ptr<bcf1_t>& record,
                          bcf1_t_plus& ans) {
    range rng(record);
    assert(rng.overlaps(site.pos));

    ans.p = record;

    ans.is_ref = is_gvcf_ref_record(record.get());

    if(bcf_get_genotypes(hdr, record.get(), &ans.gt.v, &ans.gt.capacity) != 2*record->n_sample || !ans.gt.v) {
        return Status::Failure("genotyper::preprocess_record: unexpected result from bcf_get_genotypes");
    }

    ans.allele_mapping.assign(record->n_allele, -1);
    ans.allele_mapping[0] = 0;
    ans.deletion_allele.assign(record->n_allele, false);
    ans.has_lost_allele = false;

    // map the bcf1_t alt alleles according to unification
    // checking for valid dna regex match
    string ref_al(record->d.allele[0]);
    for (int i = 1; i < record->n_allele; i++) {
        string al(record->d.allele[i]);
        if (regex_match(al, regex_dna)) {
            auto p = site.unification.find(allele(rng, al));
            if (p != site.unification.end()) {
                ans.allele_mapping[i] = p->second;
            } else {
                ans.has_lost_allele = true;
            }
        }
        if (al.size() < rng.size() && rng.size() == ref_al.size()) {
            ans.deletion_allele[i] = is_deletion(ref_al, al);
        }
    }

    return Status::OK();
}

///////////////////////////////////////////////////////////////////////////////
// Helpers for lifting over FORMAT fields from input to output VCF records
// e.g. GQ, AD, SB, etc.
///////////////////////////////////////////////////////////////////////////////

class IFormatFieldHelper {

public:
    // basic information for the retained field
    const retained_format_field& field_info;


    // Expected number of samples in output bcf record
    const int n_samples;

    // Expected number of values per sample in the **output** (ie unified site)
    // Note this may differ from the count of input for RetainedFieldNumber::ALLELES
    // and RetainedFieldNumber::ALT since the number of alleles may differ in input and
    // output
    const int count;

    IFormatFieldHelper(const retained_format_field& field_info_, int n_samples_, int count_) : field_info(field_info_), n_samples(n_samples_), count(count_) {}

    IFormatFieldHelper() = default;

    virtual Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header, bcf1_t* record,
                                   const map<int, int>& sample_mapping, const vector<int>& allele_mapping,
                                   const vector<string>& field_names, int n_val_per_sample) = 0;

    virtual Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header, bcf1_t* record,
                                   const map<int, int>& sample_mapping, const vector<int>& allele_mapping) = 0;

    virtual Status censor(int sample, bool half_call) = 0;

    virtual Status update_record_format(const bcf_hdr_t* hdr, bcf1_t* record) = 0;

    virtual ~IFormatFieldHelper() = default;
};

template <class T>
class FormatFieldHelper : public IFormatFieldHelper {
protected:
    // A vector of length n_samples * count, where each element
    // is a vector consisting of format field value from
    // record(s) related to the given sample, at a specified
    // position (for fields with more than 1 value per sample).
    // The second layer is a vector containing up to n_record values
    // whereby 0 or 1 value is pushed to it for every record processed by
    // add_record_data
    vector<vector<T>> format_v;

    // The FORMAT fields of some samples may need to be censored (emitted
    // as missing) under certain circumstances where they might otherwise
    // be unreliable/misleading. In some cases we have a flag to censor
    // only fields discussing the reference allele (for "half-calls")
    set<pair<int,bool>> censored_samples;

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
    htsvecbox<int32_t> iv_;
    htsvecbox<float> fv_;
    int bcf_get_format_wrapper(const bcf_hdr_t* dataset_header, bcf1_t* record, const char* field_name, int32_t** v) {
        int ans = bcf_get_format_int32(dataset_header, record, field_name, &iv_.v, &iv_.capacity);
        *v = iv_.v;
        return ans;
    }
    int bcf_get_format_wrapper(const bcf_hdr_t* dataset_header, bcf1_t* record, const char* field_name, float** v) {
        int ans = bcf_get_format_float(dataset_header, record, field_name, &fv_.v, &fv_.capacity);
        *v = fv_.v;
        return ans;
    }

    Status get_missing_value(T* val) {
        switch (field_info.type) {
            case RetainedFieldType::INT:
                *val = bcf_int32_missing;
                break;
            case RetainedFieldType::FLOAT:
                // Union construct to side-step compiler warnings about type checking
                // For reference, refer to bcf_float_set function in htslib/vcf.h
                union {uint32_t i; float f; } u;
                u.i = bcf_float_missing;
                *val = u.f;
                break;
            default:
                return Status::Invalid("genotyper: encountered unknown field_info.type");
        }
        return Status::OK();
    }

    Status get_default_value(T* val) {
        if (field_info.default_type == DefaultValueFiller::ZERO) {
            *val = 0;
        } else if (field_info.default_type == DefaultValueFiller::MISSING) {
            Status s;
            S(get_missing_value(val));
        } else {
            return Status::Invalid("genotyper: encountered unknown default value filler type");
        }
        return Status::OK();
    }
    Status combine_format_data(vector<T>& ans) {
        Status s;
        ans.clear();

        // Templatized default value
        T default_value;
        S(get_default_value(&default_value));

        for (auto& v : format_v) {
            if (v.empty()) {
                v.push_back(default_value);
            }
        }

        assert(format_v.size() == n_samples * count);

        // Combine values using the combine_f function given
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

    virtual Status perform_censor(vector<T>& values) {
        Status s;
        if (!censored_samples.empty()) {
            T missing_value;
            S(get_missing_value(&missing_value));
            for (auto& cs : censored_samples) {
                for (int j = 0; j < count; j++) {
                    values[cs.first*count+j] = missing_value;
                }
            }
        }
        return Status::OK();
    }

public:

    FormatFieldHelper(const retained_format_field& field_info_, int n_samples_, int count_) : IFormatFieldHelper(field_info_, n_samples_, count_) {

        switch (field_info.combi_method) {
            case FieldCombinationMethod::MIN:
                combine_f = FormatFieldHelper::min_element_wrapper;
                break;
            case FieldCombinationMethod::MAX:
                combine_f = FormatFieldHelper::max_element_wrapper;
                break;
        }

        format_v.resize(n_samples_ * count_);
    }

    virtual ~FormatFieldHelper() = default;

    // Wrapper with default values populated for
    // field_names and n_val_per_sample
    Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header,
                           bcf1_t* record, const map<int, int>& sample_mapping,
                           const vector<int>& allele_mapping) {

        return add_record_data(dataset, dataset_header, record, sample_mapping, allele_mapping, {}, -1);
    }

    Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header,
                           bcf1_t* record, const map<int, int>& sample_mapping,
                           const vector<int>& allele_mapping,
                           const vector<string>& field_names, int n_val_per_sample) {

        bool found = false;
        if (n_val_per_sample < 0) {
            n_val_per_sample = expected_n_val_per_sample(record);
        }

        const vector<string> * names_to_search = &field_names;
        if (names_to_search->empty()) {
            names_to_search = &(field_info.orig_names);
        }

        for (auto& field_name : *names_to_search) {
            if (found) break;
            T *v = nullptr;

            // rv is the number of values written
            int rv = FormatFieldHelper::bcf_get_format_wrapper(dataset_header, record, field_name.c_str(), &v);

            // raise error if there's a failed get due to type mismatch
            if (rv == -2) {
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
                        assert(in_ind < rv);
                        format_v[out_ind].push_back(v[in_ind]);
                    } // close for j loop
                } // close for i loop
            } // close rv >= 0
        }

        return found ? Status::OK() : Status::NotFound();
    }

    // Mark the (output) sample as censored.
    Status censor(int sample, bool half_call) {
        if (sample < 0 || sample >= n_samples) return Status::Invalid("genotyper::FormatFieldHelper::censor");
        censored_samples.insert(make_pair(sample, half_call));
        return Status::OK();
    }

    Status update_record_format(const bcf_hdr_t* hdr, bcf1_t* record) {
        Status s;
        vector<T> ans;
        S(combine_format_data(ans));
        assert(ans.size() == n_samples*count);
        S(perform_censor(ans));

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

// Special-case logic for the allele depth (AD) field
class ADFieldHelper : public FormatFieldHelper<int32_t> {
public:
    ADFieldHelper(const retained_format_field& field_info_, int n_samples_, int count_)
        : FormatFieldHelper<int32_t>(field_info_, n_samples_, count_) {
        assert(field_info.name == "AD");
    }

    Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header,
                           bcf1_t* record, const map<int, int>& sample_mapping,
                           const vector<int>& allele_mapping,
                           const vector<string>& field_names, int n_val_per_sample) {
        Status s = FormatFieldHelper<int32_t>::add_record_data(dataset, dataset_header, record, sample_mapping, allele_mapping, field_names, n_val_per_sample);

        if (s == StatusCode::NOT_FOUND) {
            // Record has no AD field, usually meaning it's a reference confidence record
            // (though there are exceptions, e.g. gVCF test case DP0_noAD).
            // use MIN_DP/DP as the reference allele depth
            s = FormatFieldHelper<int32_t>::add_record_data(dataset, dataset_header, record, sample_mapping, allele_mapping, {"MIN_DP", "DP"}, 1);
        }

        return s;
    }

protected:
    Status perform_censor(vector<int32_t>& values) override {
        Status s;
        if (!censored_samples.empty()) {
            int32_t missing_value;
            S(get_missing_value(&missing_value));
            for (auto& cs : censored_samples) {
                // if a half-call (cs.second), censor only the reference allele
                // depth, which can be misleading when there are overlapping or
                // unphased records.
                for (int j = 0; j < (cs.second ? 1 : count); j++) {
                    values[cs.first*count+j] = missing_value;
                }
            }
        }
        return Status::OK();
    }
};

Status setup_format_helpers(vector<unique_ptr<IFormatFieldHelper>>& format_helpers,
                            const vector<retained_format_field>& liftover_fields,
                            const unified_site& site,
                            const vector<string>& samples) {
    for (const auto& format_field_info : liftover_fields) {
        int count = -1;
        if (format_field_info.number == RetainedFieldNumber::BASIC) {
            count = format_field_info.count;
        } else if (format_field_info.number == RetainedFieldNumber::ALT) {
            // site.alleles.size() gives # alleles incl. REF
            count = (site.alleles.size() - 1);
        } else if (format_field_info.number == RetainedFieldNumber::ALLELES) {
            count = (site.alleles.size());
        }

        if (count < 0) {
            return Status::Failure("setup_format_helpers: failed to identify count for format field");
        }

        if (format_field_info.name == "AD" && format_field_info.number == RetainedFieldNumber::ALLELES) {
            format_helpers.push_back(unique_ptr<IFormatFieldHelper>(new ADFieldHelper(format_field_info, samples.size(), count)));
        } else switch (format_field_info.type) {
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

    return Status::OK();
}

Status update_format_fields(const string& dataset, const bcf_hdr_t* dataset_header, const map<int,int>& sample_mapping,
                            const unified_site& site, vector<unique_ptr<IFormatFieldHelper>>& format_helpers,
                            const vector<shared_ptr<bcf1_t_plus>>& all_records,
                            const vector<shared_ptr<bcf1_t_plus>>& variant_records) {
    Status s;

    // Update format helpers
    for (auto& format_helper : format_helpers) {
        const vector<shared_ptr<bcf1_t_plus>> *records_to_use = nullptr;
        if (format_helper->field_info.ignore_non_variants && !variant_records.empty()) {
            // Only care about variant records, loop through variant_records
            records_to_use = &variant_records;
        } else {
            // Look through all records (variant and non_variant)
            records_to_use = &all_records;
        }

        for (const auto& record : *records_to_use) {
            s = format_helper->add_record_data(dataset, dataset_header, record->p.get(),
                                               sample_mapping, record->allele_mapping);
            if (s.bad() && s != StatusCode::NOT_FOUND) {
                return s;
            }
        }
    }
    return Status::OK();
}

// Helper class for keeping track of the per-allele depth of coverage info in
// a bcf1_t record. There are a couple different cases to handle, depending on
// whether we're looking at a gVCF reference confidence record or a "regular"
// VCF record.
class AlleleDepthHelper {
    const genotyper_config& cfg_;
    size_t n_sample_ = 0, n_allele_ = 0;
    bool is_g_ = false;
    htsvecbox<int32_t> v_;

public:

    // The AlleleDepthHelper is constructed into an undefined state. Load()
    // must be invoked, successfully, before it can be used.
    AlleleDepthHelper(const genotyper_config& cfg)
        : cfg_(cfg)
        {}

    // The helper can be reused for multiple records by calling Load()
    // repeatedly. This will be slightly more efficient than using a new
    // helper for each record.
    Status Load(const string& dataset, const bcf_hdr_t* dataset_header, bcf1_t* record) {
        n_sample_ = record->n_sample;
        n_allele_ = record->n_allele;
        // is this a gVCF reference confidence record?
        is_g_ = is_gvcf_ref_record(record);

        if (is_g_) {
            // if so, look for the MIN_DP FORMAT field (or equivalent)
            int nv = bcf_get_format_int32(dataset_header, record, cfg_.ref_dp_format.c_str(),
                                          &v_.v, &v_.capacity);

            if (nv != record->n_sample) {
                ostringstream errmsg;
                errmsg << dataset << " " << range(record).str() << " (" << cfg_.ref_dp_format << ")";
                return Status::Invalid("genotyper: gVCF reference depth FORMAT field is missing or malformed", errmsg.str());
            }
        } else {
            // this is a regular VCF record, so look for the AD FORMAT field (or equivalent)
            int nv = bcf_get_format_int32(dataset_header, record, cfg_.allele_dp_format.c_str(),
                                          &v_.v, &v_.capacity);

            if (nv == -1 || nv == -3) {
                // We allow the AD field to not exist mainly as a (poor)
                // workaround for some of our test case gVCFs not having it...

                if (nv == -3) {
                    // AD is declared in the header, but not present in the
                    // FORMAT fields for this record. We'll tolerate this for
                    // an unusual observed class of variant records which have
                    // INFO DP=0 (gVCF test case DP0_noAD.yml)

                    nv = bcf_get_info_int32(dataset_header, record, "DP", &v_.v, &v_.capacity);
                    if (nv != 1 || v_[0] != 0) {
                        ostringstream errmsg;
                        errmsg << dataset << " " << range(record).str() << " (" << cfg_.allele_dp_format << ")";
                        return Status::Invalid("genotyper: VCF allele depth FORMAT field is missing", errmsg.str());
                    }
                }

                nv = record->n_sample * record->n_allele;
                if (v_.capacity < nv) {
                    v_.v = (int32_t*) realloc(v_.v, nv*sizeof(int32_t));
                    v_.capacity = nv;
                }
                memset(v_.v, 0, nv*sizeof(int32_t));
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

    bool is_gvcf_ref() { return is_g_; }

    // get depth for sample i & allele j
    unsigned get(unsigned sample, unsigned allele) {
        if (sample >= n_sample_ || allele >= n_allele_) return 0;
        if (is_g_) {
            // the MIN_DP array has just one integer per sample
            if (allele != 0) return 0;
            return v_[sample];
        } else {
            // the AD array has one integer per allele per sample
            return v_[sample*n_allele_+allele];
        }
    }
};

/// A helper function to update min_ref_depth based on several reference
/// confidence records. min_ref_depth[j] is the minimum depth of reference
/// coverage seen for sample j across the reference confidence records, and
/// should be initialized to -1 before any reference confidence records are
/// seen.
static Status update_min_ref_depth(const string& dataset, const bcf_hdr_t* dataset_header,
                                   int bcf_nsamples, const map<int,int>& sample_mapping,
                                   const vector<shared_ptr<bcf1_t_plus>>& ref_records,
                                   AlleleDepthHelper& depth,
                                   vector<int>& min_ref_depth) {
    Status s;
    for (auto& ref_record : ref_records) {
        S(depth.Load(dataset, dataset_header, ref_record->p.get()));

        for (int i=0; i<bcf_nsamples; i++) {
            int mapped_sample = sample_mapping.at(i);
            assert(mapped_sample < min_ref_depth.size());

            if (min_ref_depth[mapped_sample] < 0) {
                min_ref_depth[mapped_sample] = depth.get(i, 0);
            } else {
                min_ref_depth[mapped_sample] = min(min_ref_depth[mapped_sample],
                                                   (int) depth.get(i, 0));
            }
        }
    }
    return Status::OK();
}


///////////////////////////////////////////////////////////////////////////////
// Genotyper core
///////////////////////////////////////////////////////////////////////////////

/// Revise genotypes which are initially called with lost alleles. Frequently these
/// are low quality and "round down" to homozygous ref. Mutates record.
Status revise_genotypes(const genotyper_config& cfg, const unified_site& us, const map<int, int>& sample_mapping,
                        const bcf_hdr_t* hdr, bcf1_t_plus& vr) {
    assert(!vr.is_ref);
    if (!vr.has_lost_allele) {
        // below would be a no-op anyway
        return Status::OK();
    }

    // start by replacing the record with a duplicate, since it may not be safe to
    // mutate the "original"
    auto record = shared_ptr<bcf1_t>(bcf_dup(vr.p.get()), &bcf_destroy);
    vr.p = record;
    if (bcf_unpack(record.get(), BCF_UN_ALL)) return Status::Failure("genotyper::prepare_dataset_records bcf_unpack");
    unsigned nGT = diploid::genotypes(record->n_allele);
    range rng(record);

    // extract input genotype likelihoods and GQ
    vector<double> gll;
    Status s = diploid::bcf_get_genotype_log_likelihoods(hdr, record.get(), gll);
    if (!s.ok()) {
        return Status::Failure("genotyper::revise_genotypes: couldn't find genotype likelihoods in gVCF record", s.str());
    }
    assert(gll.size() == nGT*record->n_sample);
    htsvecbox<int32_t> gq;
    if(bcf_get_format_int32(hdr, record.get(), "GQ", &gq.v, &gq.capacity) != record->n_sample || !gq.v) {
        return Status::Failure("genotyper::revise_genotypes: unexpected result from bcf_get_format_int32 GQ");
    }
    assert(gq.capacity >= record->n_sample);

    // construct "prior" over input ALT alleles which penalizes lost ones (otherwise flat)
    const float lost_log_prior = log(std::max(us.lost_allele_frequency, cfg.min_assumed_allele_frequency));
    vector<double> gt_log_prior(diploid::genotypes(record->n_allele), 0.0);
    for (int i = 0; i < record->n_allele; i++) {
        if (vr.allele_mapping[i] == -1) {
            assert(i > 0);
            gt_log_prior[i] = lost_log_prior;
        }
    }

    // proceed through designated samples
    for (const auto& sample : sample_mapping) {
        assert(sample.first < record->n_sample);
        // add "priors" to genotype likelihoods; keep track of MAP and 2nd (silver)
        double* sample_gll = gll.data() + sample.first*nGT;
        double map_gll = log(0), silver_gll = log(0);
        int map_gt = -1;
        for (int g = 0; g < nGT; g++) {
            const auto alleles = diploid::gt_alleles(g);
            // Use the smaller of the priors on the two alleles and not their product.
            // If we view this as "penalizing" the likelihoods of genotypes which include
            // lost alleles, one such penalty is sufficient.
            auto g_ll = sample_gll[g] + std::min(gt_log_prior[alleles.first], gt_log_prior[alleles.second]);
            if (g_ll > map_gll) {
                silver_gll = map_gll;
                map_gll = g_ll;
                map_gt = g;
            } else if (g_ll > silver_gll) {
                silver_gll = g_ll;
            }
        }
        assert(map_gt >= 0 && map_gt < nGT);
        assert(map_gll >= silver_gll);
        assert(silver_gll > log(0));

        // record MAP genotype and recalculate GQ
        const auto revised_alleles = diploid::gt_alleles(map_gt);
        vr.gt.v[sample.first*2] = bcf_gt_unphased(revised_alleles.first);
        vr.gt.v[sample.first*2+1] = bcf_gt_unphased(revised_alleles.second);
        gq.v[sample.first] = std::min(99, (int) round(10.0*(map_gll - silver_gll)/log(10.0)));
    }

    // write GT and GQ back into record
    if (bcf_update_format_int32(hdr, record.get(), "GQ", gq.v, record->n_sample)) {
        return Status::Failure("genotyper::revise_genotypes: bcf_update_format_int32 GQ failed");
    }
    if (bcf_update_genotypes(hdr, record.get(), vr.gt.v, 2*record->n_sample)) {
        return Status::Failure("genotyper::revise_genotypes: bcf_update_genotypes failed");
    }

    return Status::OK();
}

/// Given a unified site and the set of gVCF records overlapping it in some
/// dataset, check that they span the site, preprocess them, and separate
/// the reference and variant records.
///
/// Ideally and often there's either zero or one variant records -- zero if
/// the dataset exhibits no variation at the site, and one if it does.
/// Unfortunately, there are a number of circumstances under which some
/// variant callers produce multiple overlapping records.
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
///      variant_records filled in, min_ref_depth updated accordingly, rnc = N_A
///
///
/// FIXME: detect & complain if the reference confidence records actually overlap the
///        variant records
Status prepare_dataset_records(const genotyper_config& cfg, const unified_site& site,
                               const string& dataset, const bcf_hdr_t* hdr, int bcf_nsamples,
                               const map<int, int>& sample_mapping,
                               const vector<shared_ptr<bcf1_t>>& records,
                               AlleleDepthHelper& depth,
                               NoCallReason& rnc,
                               vector<int>& min_ref_depth,
                               vector<shared_ptr<bcf1_t_plus>>& all_records,
                               vector<shared_ptr<bcf1_t_plus>>& variant_records) {
    // initialize outputs
    rnc = NoCallReason::MissingData;
    all_records.clear();
    variant_records.clear();

    Status s;

    // collect the ranges covered by the records
    vector<range> record_rngs;
    record_rngs.reserve(records.size());
    for (const auto& record: records) {
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

    vector<shared_ptr<bcf1_t_plus>> ref_records;
    for (const auto& record : records) {
        auto rp = make_shared<bcf1_t_plus>();
        S(preprocess_record(site, hdr, record, *rp));
        if (rp->is_ref) {
            ref_records.push_back(rp);
        } else {
            if (cfg.revise_genotypes) {
                S(revise_genotypes(cfg, site, sample_mapping, hdr, *rp));
            }
            variant_records.push_back(rp);
        }
        all_records.push_back(rp);
    }

    // compute min_ref_depth across the reference confidence records
    S(update_min_ref_depth(dataset, hdr, bcf_nsamples, sample_mapping,
                           ref_records, depth, min_ref_depth));

    // Success...
    rnc = NoCallReason::N_A;
    return Status::OK();
}

/// Based on the cluster of variant records and min_ref_depth produced by
/// prepare_dataset_records, fill genotypes for this dataset's samples with
/// appropriate calls (currently by translation of the input hard-calls).
/// Updates genotypes and may modify min_ref_depth.
///
/// FIXME: not coded to deal with multi-sample gVCFs properly.
static Status translate_genotypes(const genotyper_config& cfg, const unified_site& site,
                                  const string& dataset, const bcf_hdr_t* dataset_header,
                                  int bcf_nsamples, const map<int,int>& sample_mapping,
                                  const vector<shared_ptr<bcf1_t_plus>>& variant_records,
                                  AlleleDepthHelper& depth,
                                  vector<int>& min_ref_depth,
                                  vector<one_call>& genotypes) {
    assert(genotypes.size() == 2*min_ref_depth.size());
    assert(!site.monoallelic);
    Status s;

    // Scan the variant records to pull out those with 0/0 genotype calls
    // from those actually presenting variation
    vector<shared_ptr<bcf1_t_plus>> records_00, records_non00;
    for (const auto& a_record : variant_records) {
        assert(!a_record->is_ref);
        bool non00 = false;
        for (int i = 0; i < 2*bcf_nsamples; i++) {
            assert(a_record->gt.capacity > i);
            if (bcf_gt_is_missing(a_record->gt[i]) || bcf_gt_allele(a_record->gt[i]) != 0) {
                non00 = true;
            }
        }
        if (non00) {
            records_non00.push_back(a_record);
        } else {
            records_00.push_back(a_record);
        }
    }

    // update min_ref_depth with found 0/0 records
    S(update_min_ref_depth(dataset, dataset_header, bcf_nsamples, sample_mapping,
                           records_00, depth, min_ref_depth));

    bcf1_t_plus* record = nullptr;
    bool half_call = false;

    if (records_non00.size() == 0) {
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
    } else if (records_non00.size() == 1)  {
        // simple common case: one variant record overlapping the unified site
        record = records_non00[0].get();
    } else {
        // complex situation: multiple non-0/0 records overlapping the unified site.
        //
        // If the records don't all share at least one position in common (i.e. their
        // ALT alleles aren't mutually exclusive on one chromosome), punt with
        // UnphasedVariants. We'll improve this in the future.
        //
        // If at least one record is a heterozygous 0/X call where X is a known
        // allele in the unified site, and none of the records call >1 ALT allele,
        // generate a half-call from the highest-quality such record.
        // This at least recovers some of the information when the GVCF has two
        // overlapping 0/X records for one sample (we'd rather it present one record
        // heterozygous for two ALTs)
        //
        // Otherwise: punt with OverlappingVariants

        half_call = true;
        range intersection(records_non00[0]->p);
        for (auto& a_record : records_non00) {
            range record_rng(a_record->p);
            assert(record_rng.rid == intersection.rid);
            intersection.beg = max(record_rng.beg, intersection.beg);
            intersection.end = min(record_rng.end, intersection.end);

            for (int i = 0; half_call && i < bcf_nsamples; i++) {
                assert(a_record->gt.capacity > 2*i);
                if (!bcf_gt_is_missing(a_record->gt[2*i]) && bcf_gt_allele(a_record->gt[2*i]) != 0) {
                    half_call = false;
                } else if (!bcf_gt_is_missing(a_record->gt[2*i+1])) {
                    auto al = bcf_gt_allele(a_record->gt[2*i+1]);
                    if (al > 0 && a_record->allele_mapping[al] > 0 &&
                        (!record || record->p->qual < a_record->p->qual)) {
                        record = a_record.get();
                    }
                }
            }
        }

        if (intersection.beg >= intersection.end) {
            for (int i = 0; i < bcf_nsamples; i++) {
                genotypes[sample_mapping.at(i)*2].RNC =
                    genotypes[sample_mapping.at(i)*2+1].RNC =
                        NoCallReason::UnphasedVariants;
            }
            return Status::OK();
        }

        if (!record || !half_call) {
            for (int i = 0; i < bcf_nsamples; i++) {
                genotypes[sample_mapping.at(i)*2].RNC =
                    genotypes[sample_mapping.at(i)*2+1].RNC =
                        NoCallReason::OverlappingVariants;
            }
            return Status::OK();
        }
    }

    // Now, translating genotypes from one variant BCF record.
    assert(record != nullptr);

    // get the genotype calls
    htsvecbox<int> gt;
    int nGT = bcf_get_genotypes(dataset_header, record->p.get(), &gt.v, &gt.capacity);
    int n_bcf_samples = bcf_hdr_nsamples(dataset_header);
    if (!gt.v || nGT != 2*n_bcf_samples) return Status::Failure("genotyper::translate_genotypes bcf_get_genotypes");
    assert(record->p->n_sample == bcf_hdr_nsamples(dataset_header));

    S(depth.Load(dataset, dataset_header, record->p.get()));

    // for each shared sample, record the genotype call.
    for (const auto& ij : sample_mapping) {
        assert(2*ij.first < nGT);
        assert(ij.second < min_ref_depth.size());

        // TODO: are depth and allele_mapping checks inside-out????
        #define fill_allele(ofs)                                                  \
            if (gt[2*ij.first+ofs] != bcf_int32_vector_end &&                     \
                !bcf_gt_is_missing(gt[2*ij.first+(ofs)])) {                       \
                auto al = bcf_gt_allele(gt[2*ij.first+(ofs)]);                    \
                assert(al >= 0 && al < record->p->n_allele);                      \
                int rd = min_ref_depth[ij.second];                                \
                if (depth.get(ij.first, al) >= cfg.required_dp                    \
                    && (rd < 0 || rd >= cfg.required_dp)) {                       \
                    if (record->allele_mapping[al] >= 0) {                        \
                        genotypes[2*ij.second+(ofs)] =                            \
                            one_call(bcf_gt_unphased(record->allele_mapping[al]), \
                                     NoCallReason::N_A);                          \
                    } else {                                                      \
                        genotypes[2*ij.second+(ofs)].RNC =                        \
                            record->deletion_allele[al]                           \
                                ? NoCallReason::LostDeletion                      \
                                : NoCallReason::LostAllele;                       \
                    }                                                             \
                } else {                                                          \
                    genotypes[2*ij.second+(ofs)].RNC =                            \
                        NoCallReason::InsufficientDepth;                          \
                }                                                                 \
            }

        if (half_call) {
            genotypes[2*ij.second].RNC = NoCallReason::OverlappingVariants;
        } else {
            fill_allele(0)
        }
        fill_allele(1)
    }

    return Status::OK();
}

/// streamlined version of translate_genotypes for monoallelic sites
/// FIXME: not coded to deal with multi-sample gVCFs properly.
static Status translate_monoallelic(const genotyper_config& cfg, const unified_site& site,
                                    const string& dataset, const bcf_hdr_t* dataset_header,
                                    int bcf_nsamples, const map<int,int>& sample_mapping,
                                    const vector<shared_ptr<bcf1_t_plus>>& variant_records,
                                    AlleleDepthHelper& depth,
                                    vector<int>& min_ref_depth,
                                    vector<one_call>& genotypes) {
    assert(genotypes.size() == 2*min_ref_depth.size());
    assert(site.monoallelic);
    assert(site.alleles.size() == 2);
    Status s;
    bcf1_t_plus* record = nullptr;

    // scan the (potentially) multiple overlapping records to find the one with
    // the desired allele
    for (auto& a_record : variant_records) {
        assert(!a_record->is_ref);
        for (unsigned al = 1; al < a_record->allele_mapping.size(); al++) {
            if (a_record->allele_mapping[al] == 1) {
                if (record == nullptr) {
                    record = a_record.get();
                } else {
                    // uh, overlapping records with the same allele!?
                    for (int i = 0; i < bcf_nsamples; i++) {
                        genotypes[sample_mapping.at(i)*2].RNC =
                            genotypes[sample_mapping.at(i)*2+1].RNC =
                                NoCallReason::OverlappingVariants;
                    }
                    return Status::OK();
                }
            }
        }
    }

    if (record == nullptr) {
        // we have nothing to say here
        for (int i = 0; i < bcf_nsamples; i++) {
            genotypes[sample_mapping.at(i)*2].RNC =
                genotypes[sample_mapping.at(i)*2+1].RNC =
                    NoCallReason::MonoallelicSite;
        }
        return Status::OK();
    }

    // Now, translating genotypes from one variant BCF record.

    // get the genotype calls
    htsvecbox<int> gt;
    int nGT = bcf_get_genotypes(dataset_header, record->p.get(), &gt.v, &gt.capacity);
    int n_bcf_samples = bcf_hdr_nsamples(dataset_header);
    if (!gt.v || nGT != 2*n_bcf_samples) return Status::Failure("genotyper::translate_genotypes bcf_get_genotypes");
    assert(record->p->n_sample == bcf_hdr_nsamples(dataset_header));

    S(depth.Load(dataset, dataset_header, record->p.get()));

    // for each shared sample, record the genotype call.
    for (const auto& ij : sample_mapping) {
        assert(2*ij.first < nGT);
        assert(ij.second < min_ref_depth.size());

        #define fill_monoallelic(ofs)                                             \
            if (gt[2*ij.first+ofs] != bcf_int32_vector_end &&                     \
                !bcf_gt_is_missing(gt[2*ij.first+(ofs)])) {                       \
                auto al = bcf_gt_allele(gt[2*ij.first+(ofs)]);                    \
                assert(al >= 0 && al < record->p->n_allele);                      \
                if (record->allele_mapping[al] > 0) {                             \
                    if (depth.get(ij.first, al) >= cfg.required_dp) {             \
                        genotypes[2*ij.second+(ofs)] =                            \
                            one_call(bcf_gt_unphased(record->allele_mapping[al]), \
                                     NoCallReason::N_A);                          \
                    } else {                                                      \
                        genotypes[2*ij.second+(ofs)].RNC =                        \
                            NoCallReason::InsufficientDepth;                      \
                    }                                                             \
                } else {                                                          \
                    genotypes[2*ij.second+(ofs)].RNC =                            \
                        NoCallReason::MonoallelicSite;                            \
                }                                                                 \
            }

        fill_monoallelic(0)
        fill_monoallelic(1)
    }

    return Status::OK();
}

Status genotype_site(const genotyper_config& cfg, MetadataCache& cache, BCFData& data, const unified_site& site,
                     const std::string& sampleset, const vector<string>& samples,
                     const bcf_hdr_t* hdr, shared_ptr<bcf1_t>& ans,
                     bool residualsFlag, shared_ptr<string> &residual_rec,
                     atomic<bool>* ext_abort) {
    Status s;

    // Initialize a vector for the unified genotype calls for each sample,
    // starting with everything missing. We'll then loop through BCF records
    // overlapping this site and fill in the genotypes as we encounter them.
    vector<one_call> genotypes(2*samples.size());

    // Setup format field helpers
    vector<unique_ptr<IFormatFieldHelper>> format_helpers;
    S(setup_format_helpers(format_helpers, cfg.liftover_fields, site, samples));

    shared_ptr<const set<string>> samples2, datasets;
    vector<unique_ptr<RangeBCFIterator>> iterators;
    S(data.sampleset_range(cache, sampleset, site.pos, nullptr,
                           samples2, datasets, iterators));
    assert(samples.size() == samples2->size());

    AlleleDepthHelper adh(cfg);
    vector<DatasetResidual> lost_calls_info;

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

        // pre-process the records
        vector<int> min_ref_depth(samples.size(), -1);
        vector<shared_ptr<bcf1_t_plus>> all_records, variant_records;
        NoCallReason rnc = NoCallReason::MissingData;
        S(prepare_dataset_records(cfg, site, dataset, dataset_header.get(), bcf_nsamples,
                                  sample_mapping, records, adh, rnc, min_ref_depth,
                                  all_records, variant_records));

        if (rnc != NoCallReason::N_A) {
            // no call for the samples in this dataset (several possible
            // reasons)
            for (const auto& p : sample_mapping) {
                genotypes[p.second*2].RNC =
                    genotypes[p.second*2+1].RNC = rnc;
            }
        } else if (!site.monoallelic) {
            // make genotype calls for the samples in this dataset
            S(translate_genotypes(cfg, site, dataset, dataset_header.get(), bcf_nsamples,
                                  sample_mapping, variant_records, adh, min_ref_depth,
                                  genotypes));
        } else {
            S(translate_monoallelic(cfg, site, dataset, dataset_header.get(), bcf_nsamples,
                                    sample_mapping, variant_records, adh, min_ref_depth,
                                    genotypes));
        }

        // Update FORMAT fields for this dataset.
        S(update_format_fields(dataset, dataset_header.get(), sample_mapping, site, format_helpers,
                               all_records, variant_records));
        // But if rnc = MissingData, PartialData, UnphasedVariants, or OverlappingVariants, then
        // we must censor the FORMAT fields as potentially unreliable/misleading.
        for (const auto& p : sample_mapping) {
            auto rnc1 = genotypes[p.second*2].RNC;
            auto rnc2 = genotypes[p.second*2+1].RNC;

            if (rnc1 == NoCallReason::MissingData || rnc1 == NoCallReason::PartialData) {
                assert(rnc1 == rnc2);
                for (const auto& fh : format_helpers) {
                    S(fh->censor(p.second, false));
                }
            } else if (site.monoallelic || rnc1 == NoCallReason::UnphasedVariants || rnc1 == NoCallReason::OverlappingVariants) {
                const bool half_call = site.monoallelic || ((genotypes[p.second*2].RNC == NoCallReason::N_A) != (genotypes[p.second*2+1].RNC == NoCallReason::N_A));

                for (const auto& fh : format_helpers) {
                    if (fh->field_info.name != "DP") {
                        S(fh->censor(p.second, half_call));
                    }
                }
            }
        }

        // Handle residuals
        if (residualsFlag) {
            const set<NoCallReason> non_residual_RNCs = { NoCallReason::N_A, NoCallReason::MissingData,
                                                          NoCallReason::PartialData, NoCallReason::InsufficientDepth };

            bool any_lost_calls = false;
            for (int i = 0; i < bcf_nsamples; i++) {
                if (non_residual_RNCs.find(genotypes[sample_mapping.at(i)*2].RNC) == non_residual_RNCs.end() ||
                    non_residual_RNCs.find(genotypes[sample_mapping.at(i)*2 + 1].RNC) == non_residual_RNCs.end()) {
                    any_lost_calls = true;
                    break;
                }
            }

            if (any_lost_calls) {
                // missing call, keep it in memory
                DatasetResidual dsr;
                dsr.name = dataset;
                dsr.header = dataset_header;
                dsr.records = records;
                lost_calls_info.push_back(dsr);
            }
        }
    }

    // Clean up emission order of alleles
    for(size_t i=0; i < samples.size(); i++) {
// TODO
//        if(genotypes[2*i].allele != bcf_gt_missing && genotypes[2*i+1].allele == bcf_gt_missing ||
//           genotypes[2*i].allele > genotypes[2*i+1].allele) {
        if(genotypes[2*i] > genotypes[2*i+1]) {
            swap(genotypes[2*i], genotypes[2*i+1]);
        }
    }
    // Create the destination BCF record for this site.
    ans = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    ans->rid = site.pos.rid;
    ans->pos = site.pos.beg;
    ans->rlen = site.pos.end - site.pos.beg;
    ans->qual = site.qual;

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

    // Lifted-over FORMAT fields (non-genotype based)
    for (auto& format_helper : format_helpers) {
        S(format_helper->update_record_format(hdr, ans.get()));
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
            RNC_CASE(MonoallelicSite,"M")
            default:
                assert(c.RNC == NoCallReason::MissingData);
        }
        rnc.push_back(v);
    }
    assert (gt.size() == rnc.size());
    if (bcf_update_format_string(hdr, ans.get(), "RNC", rnc.data(), rnc.size()) != 0) {
        return Status::Failure("bcf_update_format_string RNC");
    }

    if (site.monoallelic && bcf_add_filter(hdr, ans.get(), bcf_hdr_id2int(hdr, BCF_DT_ID, "MONOALLELIC")) != 1) {
        return Status::Failure("bcf_add_filter MONOALLELIC");
    }

    if (residualsFlag &&
        !lost_calls_info.empty()) {
        // Write loss record to the residuals file, useful for offline debugging.
        residual_rec = make_shared<string>();
        S(residuals_gen_record(site, hdr, ans.get(), lost_calls_info,
                               cache, samples,
                               *residual_rec));
    }

    // Overwrite the output BCF record with a duplicate. Why? This forces htslib to
    // perform some internal serialization of the data (see the static bcf1_sync
    // function in vcf.c, which we can't call directly, but is called by bcf_dup).
    // htslib would otherwise do this serialization implicitly while writing the
    // record out to a file, but by doing it explicitly here, we get to do some of the
    // work in the current worker thread rather than the single thread responsible for
    // writing out the file.
    auto ans2 = shared_ptr<bcf1_t>(bcf_dup(ans.get()), &bcf_destroy);
    ans = move(ans2);

    return Status::OK();
}

}
