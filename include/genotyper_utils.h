// Helper classes/functions for the genotyper algorithm (included only by genotyper.cc)
namespace GLnexus {

///////////////////////////////////////////////////////////////////////////////
// Helpers for lifting over FORMAT fields from input to output VCF records
// e.g. GQ, AD, SB, etc.
///////////////////////////////////////////////////////////////////////////////

class FormatFieldHelper {

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

    FormatFieldHelper(const retained_format_field& field_info_, int n_samples_, int count_) : field_info(field_info_), n_samples(n_samples_), count(count_) {}

    FormatFieldHelper() = default;

    virtual Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header, bcf1_t* record,
                                   const map<int, int>& sample_mapping, const vector<int>& allele_mapping,
                                   const int n_allele_out, const vector<string>& field_names, int n_val_per_sample) = 0;

    // Wrapper with default values populated for
    // field_names and n_val_per_sample
    virtual Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header, bcf1_t* record,
                                   const map<int, int>& sample_mapping, const vector<int>& allele_mapping,
                                   const int n_allele_out) {
        return add_record_data(dataset, dataset_header, record, sample_mapping, allele_mapping, n_allele_out, {}, -1);
    }

    virtual Status censor(int sample, bool half_call) {
        if (sample < 0 || sample >= n_samples) return Status::Invalid("genotyper::FormatFieldHelper::censor");
        censored_samples.insert(make_pair(sample, half_call));
        return Status::OK();
    }

    virtual Status update_record_format(const bcf_hdr_t* hdr, bcf1_t* record) = 0;

    virtual ~FormatFieldHelper() = default;

protected:

    // The FORMAT fields of some samples may need to be censored (emitted
    // as missing) under certain circumstances where they might otherwise
    // be unreliable/misleading. In some cases we have a flag to censor
    // only fields discussing the reference allele (for "half-calls")
    set<pair<int,bool>> censored_samples;

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
            auto n = max(record->n_allele, 2U); // accommodate GVCF ALT=. representation for non-variants
            expected_count = diploid::genotypes(n);
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
                                const vector<int>& allele_mapping,
                                const int n_allele_out) {
        int mapped_i = sample_mapping.at(unmapped_i);

        // Sample should always be mappable
        assert (mapped_i >= 0);

        switch (field_info.number){
            case RetainedFieldNumber::ALT:
            case RetainedFieldNumber::ALLELES:
            {
                // Fall through for both cases that require allele_mapping
                int mapped_j = allele_mapping.at(unmapped_j);
                if (mapped_j < 0) {
                    // Allele is not mappable (trimmed or is a gvcf record)
                    return -1;
                } else {
                    return mapped_i * count + mapped_j;
                }
            }
            case RetainedFieldNumber::GENOTYPE:
            {
                // One value per diploid genotype.
                pair<unsigned,unsigned> unmapped_alleles = diploid::gt_alleles(unmapped_j);
                int mapped_al1, mapped_al2;
                if (allele_mapping.size() > 1) {
                    mapped_al1 = allele_mapping.at(unmapped_alleles.first);
                    mapped_al2 = allele_mapping.at(unmapped_alleles.second);
                } else {
                    mapped_al1 = unmapped_alleles.first <= 1 ? unmapped_alleles.first : -1;
                    mapped_al2 = unmapped_alleles.second <= 1 ? unmapped_alleles.second : -1;
                }
                assert(mapped_al1 < n_allele_out && mapped_al2 < n_allele_out);
                if (mapped_al1 < 0 || mapped_al2 < 0) {
                    // TODO: arguably PL should be censored if the max-likelihood entry (0) does
                    // not map into the output vector. They may be open to misinterpretation
                    // otherwise.
                    return -1;
                }
                int mapped_j = diploid::alleles_gt(mapped_al1, mapped_al2);
                assert(mapped_j >= 0 && mapped_j < count);
                return mapped_i * count + mapped_j;
            }
            default:
            {
                // RetainedFieldNumber::BASIC case
                return mapped_i * count + unmapped_j;
            }
        }
    }
};

template <class T>
class NumericFormatFieldHelper : public FormatFieldHelper {
protected:
    // A vector of length n_samples * count, where each element
    // is a vector consisting of format field value from
    // record(s) related to the given sample, at a specified
    // position (for fields with more than 1 value per sample).
    // The second layer is a vector containing up to n_record values
    // whereby 0 or 1 value is pushed to it for every record processed by
    // add_record_data
    vector<vector<T>> format_v;

    // Combination function to handle combining multiple format values
    // from multiple records
    T (*combine_f) (const vector<T>&, T);

    static T max_element_wrapper(const vector<T>& v, T missing) {
        return (*max_element(v.begin(), v.end()));
    }

    static T min_element_wrapper(const vector<T>& v, T missing) {
        return (*min_element(v.begin(), v.end()));
    }

    static T missing_element_wrapper(const vector<T>& v, T missing) {
        return v.size() == 1 ? v[0] : missing;
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

    Status get_missing_value(int32_t& val) {
        assert(field_info.type == RetainedFieldType::INT);
        val = bcf_int32_missing;
        return Status::OK();
    }

    Status get_missing_value(float& val) {
        assert(field_info.type == RetainedFieldType::FLOAT);
        // Union construct to side-step compiler warnings about type checking
        // For reference, refer to bcf_float_set function in htslib/vcf.h
        union {uint32_t i; float f; } u;
        u.i = bcf_float_missing;
        val = u.f;
        return Status::OK();
    }

    Status get_missing_value(string& val) {
        assert(field_info.type == RetainedFieldType::STRING);
        val.clear();
        return Status::OK();
    }

    Status get_default_value(T& val) {
        if (field_info.default_type == DefaultValueFiller::ZERO) {
            val = 0;
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

        // Templatized missing & default values
        T missing_value, default_value;
        S(get_missing_value(missing_value));
        S(get_default_value(default_value));

        for (auto& v : format_v) {
            if (v.empty()) {
                v.push_back(default_value);
            }
        }

        assert(format_v.size() == n_samples * count);

        // Combine values using the combine_f function given
        for (auto& format_one : format_v) {
            ans.push_back(combine_f(format_one, missing_value));
        }

        return Status::OK();
    }

    virtual Status perform_censor(vector<T>& values) {
        Status s;
        if (!censored_samples.empty()) {
            T missing_value;
            S(get_missing_value(missing_value));
            for (auto& cs : censored_samples) {
                for (int j = 0; j < count; j++) {
                    values[cs.first*count+j] = missing_value;
                }
            }
        }
        return Status::OK();
    }

public:

    NumericFormatFieldHelper(const retained_format_field& field_info_, int n_samples_, int count_) : FormatFieldHelper(field_info_, n_samples_, count_) {

        switch (field_info.combi_method) {
            case FieldCombinationMethod::MIN:
                combine_f = min_element_wrapper;
                break;
            case FieldCombinationMethod::MAX:
                combine_f = max_element_wrapper;
                break;
            default:
                combine_f = missing_element_wrapper;
                break;
        }

        format_v.resize(n_samples_ * count_);
    }

    virtual ~NumericFormatFieldHelper() = default;

    Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header,
                           bcf1_t* record, const map<int, int>& sample_mapping,
                           const vector<int>& allele_mapping, const int n_allele_out,
                           const vector<string>& field_names, int n_val_per_sample) override {

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
            int rv = bcf_get_format_wrapper(dataset_header, record, field_name.c_str(), &v);

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
                    errmsg << dataset << " " << range(record).str() << " " << field_name << " vector length " << rv << ", expected " << record->n_sample * n_val_per_sample;
                    return Status::Invalid("genotyper: unexpected result when fetching record FORMAT field", errmsg.str());
                } // close rv != record->n_sample * count

                for (int i=0; i<record->n_sample; i++) {
                    for (int j=0; j<n_val_per_sample; j++) {

                        int in_ind = i * n_val_per_sample + j;
                        int out_ind = get_out_ind_of_value(i, j, sample_mapping, allele_mapping, n_allele_out);

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

    Status update_record_format(const bcf_hdr_t* hdr, bcf1_t* record) override {
        Status s;
        vector<T> ans;
        S(combine_format_data(ans));
        assert(ans.size() == n_samples*count);
        S(perform_censor(ans));

        int retval  = 0;
        switch (field_info.type) {
            case RetainedFieldType::INT:
                for (int i = 0; i < n_samples; i++) {
                    bool all_missing = true;
                    for (int j = 0; j < count; j++) {
                        if (ans[i*count+j] != bcf_int32_missing) {
                            all_missing = false;
                        }
                    }
                    if (all_missing && field_info.number != RetainedFieldNumber::BASIC) {
                        // If all entries for a sample are missing, make the VCF output field have just
                        // one missing value (.) instead of an array of missing values (eg .,.,.,.) to
                        // save space in variable-length vectors.
                        ans[i*count+1] = bcf_int32_vector_end;
                    }
                }
                retval = bcf_update_format_int32(hdr, record, field_info.name.c_str(), ans.data(), n_samples * count);
                break;
            case RetainedFieldType::FLOAT:
                retval = bcf_update_format_float(hdr, record, field_info.name.c_str(), ans.data(), n_samples * count);
                break;
            default:
                return Status::Invalid("genotyper: Unexpected RetainedFieldType when executing update_record_format.", field_info.name);
        }
        if (retval != 0) {
            return Status::Failure("genotyper: failed to update record format when executing update_record_format.", field_info.name);
        }
        return Status::OK();
    }
};

// Special-case logic for the allele depth (AD) field
class ADFieldHelper : public NumericFormatFieldHelper<int32_t> {
public:
    ADFieldHelper(const string& ref_dp_format, const retained_format_field& field_info_, int n_samples_, int count_)
        : NumericFormatFieldHelper<int32_t>(field_info_, n_samples_, count_), ref_dp_format_(ref_dp_format) {
        assert(field_info.name == "AD");
    }

    Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header,
                           bcf1_t* record, const map<int, int>& sample_mapping,
                           const vector<int>& allele_mapping, const int n_allele_out,
                           const vector<string>& field_names, int n_val_per_sample) override {
        Status s = NumericFormatFieldHelper<int32_t>::add_record_data(dataset, dataset_header, record, sample_mapping, allele_mapping, n_allele_out, field_names, n_val_per_sample);

        if (s == StatusCode::NOT_FOUND) {
            // Record has no AD field, usually meaning it's a reference confidence record
            // (though there are exceptions, e.g. gVCF test case DP0_noAD).
            // use MIN_DP/DP as the reference allele depth
            s = NumericFormatFieldHelper<int32_t>::add_record_data(dataset, dataset_header, record, sample_mapping, allele_mapping, n_allele_out, {ref_dp_format_}, 1);
        }

        return s;
    }

protected:
    const string& ref_dp_format_;

    Status perform_censor(vector<int32_t>& values) override {
        Status s;
        if (!censored_samples.empty()) {
            int32_t missing_value;
            S(get_missing_value(missing_value));
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

class StringFormatFieldHelper : public FormatFieldHelper {
protected:
    vector<vector<string>> format_v;

    Status combine_format_data(vector<string>& ans) {
        Status s;
        ans.clear();

        assert(format_v.size() == n_samples * count);

        for (auto& format_one : format_v) {
            if (format_one.size() == 1) {
                ans.push_back(format_one[0]);
            } else if (format_one.empty() || field_info.combi_method == FieldCombinationMethod::MISSING) {
                ans.push_back(".");
            } else if (field_info.combi_method == FieldCombinationMethod::SEMICOLON) {
                ostringstream oss;
                bool first = true;
                for (auto& s : format_one) {
                    if (!first) {
                        oss << ';';
                    }
                    oss << s;
                    first = false;
                }
                ans.push_back(oss.str());
            } else {
                return Status::Invalid("genotyper misconfiguration: unsupported combi_method for string format field.", field_info.name);
            }
        }

        return Status::OK();
    }

    virtual Status perform_censor(vector<string>& values) {
        Status s;
        if (!censored_samples.empty()) {
            for (auto& cs : censored_samples) {
                for (int j = 0; j < count; j++) {
                    values[cs.first*count+j] = ".";
                }
            }
        }
        return Status::OK();
    }

public:

    StringFormatFieldHelper(const retained_format_field& field_info_, int n_samples_, int count_)
        : FormatFieldHelper(field_info_, n_samples_, count_) {
            format_v.resize(n_samples_ * count_);
        }

    virtual ~StringFormatFieldHelper() = default;

    Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header,
                           bcf1_t* record, const map<int, int>& sample_mapping,
                           const vector<int>& allele_mapping, const int n_allele_out,
                           const vector<string>& field_names, int n_val_per_sample) override {

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
            char **v_raw = nullptr;
            int vsz = 0;
            std::shared_ptr<char*> v;

            int rv = bcf_get_format_string(dataset_header, record, field_name.c_str(), &v_raw, &vsz);
            if (rv >= 0) {
                v = shared_ptr<char*>(v_raw, [&v_raw](char**) { free(v_raw[0]); free(v_raw); });
            }

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
                }

                for (int i=0; i<record->n_sample; i++) {
                    for (int j=0; j<n_val_per_sample; j++) {
                        int in_ind = i * n_val_per_sample + j;
                        int out_ind = get_out_ind_of_value(i, j, sample_mapping, allele_mapping, n_allele_out);

                        if (out_ind < 0) {
                           continue;
                        }

                        assert(out_ind < format_v.size());
                        assert(in_ind < rv);
                        format_v[out_ind].push_back(v_raw[in_ind]);
                    }
                }
            }
        }

        return found ? Status::OK() : Status::NotFound();
    }

    Status update_record_format(const bcf_hdr_t* hdr, bcf1_t* record) override {
        Status s;
        vector<string> ans;
        S(combine_format_data(ans));
        assert(ans.size() == n_samples*count);
        S(perform_censor(ans));

        vector<const char*> cstrs;
        for (auto& s : ans) {
            cstrs.push_back(s.c_str());
        }
        assert(cstrs.size() == n_samples*count);
        //cstrs = {"foo","bar"};

        int retval  = bcf_update_format_string(hdr, record, field_info.name.c_str(), cstrs.data(), n_samples*count);
        if (retval != 0) {
            return Status::Failure("genotyper: failed to update record format when executing update_record_format.", field_info.name);
        }
        return Status::OK();
    }
};


// special helper to move the FILTER column of input gVCF into a FORMAT field of the output pVCF (FT)
class FilterFormatFieldHelper : public StringFormatFieldHelper {
public:
    FilterFormatFieldHelper(const retained_format_field& field_info_, int n_samples_, int count_)
        : StringFormatFieldHelper(field_info_, n_samples_, count_)
        {}

    virtual ~FilterFormatFieldHelper() = default;

    Status add_record_data(const string& dataset, const bcf_hdr_t* dataset_header,
                            bcf1_t* record, const map<int, int>& sample_mapping,
                            const vector<int>& allele_mapping, const int n_allele_out,
                            const vector<string>& field_names, int n_val_per_sample) override {
        if (n_val_per_sample < 0) {
            n_val_per_sample = expected_n_val_per_sample(record);
        }

        if (record->d.n_flt) {
            ostringstream filters;
            bool first=true;
            for (int i = 0; i < record->d.n_flt; i++) {
                string filter = dataset_header->id[BCF_DT_ID][record->d.flt[i]].key;
                if (filter.size() && filter != "PASS") {
                    // TODO: scheme to deduplicate filters seen in multiple records
                    if (!first) {
                        filters << ';';
                    }
                    filters << filter;
                    first=false;
                }
            }

            if (!first) {
                for (int i = 0; i < record->n_sample; i++) {
                    int out_ind = get_out_ind_of_value(i, 0, sample_mapping, allele_mapping, n_allele_out);
                    if (out_ind >= 0) {
                        assert(out_ind < format_v.size());
                        format_v[out_ind].push_back(filters.str());
                    }
                }
            }
        }

        return Status::OK();
    }
};


Status setup_format_helpers(vector<unique_ptr<FormatFieldHelper>>& format_helpers,
                            const genotyper_config& cfg,
                            const unified_site& site,
                            const vector<string>& samples) {
    for (const auto& format_field_info : cfg.liftover_fields) {
        int count = -1;
        if (format_field_info.number == RetainedFieldNumber::BASIC) {
            count = format_field_info.count;
        } else if (format_field_info.number == RetainedFieldNumber::ALT) {
            // site.alleles.size() gives # alleles incl. REF
            count = (site.alleles.size() - 1);
        } else if (format_field_info.number == RetainedFieldNumber::ALLELES) {
            count = (site.alleles.size());
        } else if (format_field_info.number == RetainedFieldNumber::GENOTYPE) {
            count = diploid::genotypes(site.alleles.size());
            // TODO: censor if count > 15 (5 alleles) to prevent explosion
        }

        if (count < 0) {
            return Status::Failure("setup_format_helpers: failed to identify count for format field");
        }

        if (format_field_info.name == "AD") {
            if (format_field_info.type != RetainedFieldType::INT || format_field_info.number != RetainedFieldNumber::ALLELES) {
                return Status::Invalid("genotyper misconfiguration: AD format field should have type=int, number=alleles");
            }
            format_helpers.push_back(unique_ptr<FormatFieldHelper>(new ADFieldHelper(cfg.ref_dp_format, format_field_info, samples.size(), count)));
        } else if (format_field_info.name == "FT") {
            if (format_field_info.type != RetainedFieldType::STRING || format_field_info.number != RetainedFieldNumber::BASIC || format_field_info.count != 1) {
                return Status::Invalid("genotyper misconfiguration: FT format field should have type=string, number=basic, count=1");
            }
            format_helpers.push_back(unique_ptr<FormatFieldHelper>(new FilterFormatFieldHelper(format_field_info, samples.size(), count)));
        } else switch (format_field_info.type) {
            case RetainedFieldType::INT:
            {
                format_helpers.push_back(unique_ptr<FormatFieldHelper>(new NumericFormatFieldHelper<int32_t>(format_field_info, samples.size(), count)));
                break;
            }
            case RetainedFieldType::FLOAT:
            {
                format_helpers.push_back(unique_ptr<FormatFieldHelper>(new NumericFormatFieldHelper<float>(format_field_info, samples.size(), count)));
                break;
            }
            case RetainedFieldType::STRING:
            {
                format_helpers.push_back(unique_ptr<FormatFieldHelper>(new StringFormatFieldHelper(format_field_info, samples.size(), count)));
                break;
            }
        }
    }

    return Status::OK();
}

Status update_format_fields(const string& dataset, const bcf_hdr_t* dataset_header, const map<int,int>& sample_mapping,
                            const unified_site& site, vector<unique_ptr<FormatFieldHelper>>& format_helpers,
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
                                               sample_mapping, record->allele_mapping, site.alleles.size());
            if (s.bad() && s != StatusCode::NOT_FOUND) {
                return s;
            }
        }
    }
    return Status::OK();
}

// Helper class for keeping track of the per-allele depth of coverage info in
// a bcf1_t record. There are a couple different cases to handle, depending on
// the upstream variant caller and whether we're looking at a reference confidence
// or variant record.
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
}
