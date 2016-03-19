#include "types.h"
#include <algorithm>
#include <regex>

using namespace std;

namespace GLnexus {

regex regex_dna     ("[ACGTN]+")
    , regex_id      ("[-_a-zA-Z0-9\\.]{1,100}")
    ;

// Add src alleles to dest alleles. Identical alleles alleles are merged,
// using the sum of their copy_numbers
Status merge_discovered_alleles(const discovered_alleles& src, discovered_alleles& dest) {
    for (auto& dsal : src) {
        UNPAIR(dsal,allele,ai)
        auto p = dest.find(allele);
        if (p == dest.end()) {
            dest[allele] = ai;
        } else {
            if (ai.is_ref != p->second.is_ref) {
                return Status::Invalid("allele appears as both REF and ALT", allele.dna + "@" + allele.pos.str());
            }
            p->second.copy_number += ai.copy_number;
        }
    }

    return Status::OK();
}


Status merge_loss_stats(const consolidated_loss& src, consolidated_loss& dest) {

    Status s;
    for (auto& sample_stats : src) {
        UNPAIR(sample_stats, sample, src_stats)
        auto target = dest.find(sample);
        if (target == dest.end()) {
            // sample not found in destination
            dest.insert(make_pair(sample, src_stats));
        } else {
            target->second += src_stats;
        }
    }
    return Status::OK();
}


Status range_yaml(const std::vector<std::pair<std::string,size_t> >& contigs,
                  const range& r, YAML::Emitter& yaml, bool omit_ref=false) {
    if (!omit_ref && (r.rid < 0 || r.rid >= contigs.size())) {
        return Status::NotFound("range_yaml: invalid rid", r.str());
    }
    yaml << YAML::Flow;
    yaml << YAML::BeginMap;
    if (!omit_ref) {
        yaml << YAML::Key << "ref" << YAML::Value << contigs[r.rid].first;
    }
    yaml << YAML::Key << "beg" << YAML::Value << r.beg+1;
    yaml << YAML::Key << "end" << YAML::Value << r.end;
    yaml << YAML::EndMap;
    return Status::OK();
}

Status range_of_yaml(const YAML::Node& yaml, const vector<pair<string,size_t> >& contigs,
                     range& ans, int default_rid = -1) {
    #define V(pred,msg) if (!(pred)) return Status::Invalid("range_of_yaml: " msg)

    V(yaml.IsMap(), "not a map at top level");

    int rid;
    const auto n_ref = yaml["ref"];
    if (n_ref) {
        V(n_ref.IsScalar(), "invalid 'ref' field");
        const string& ref = n_ref.Scalar();
        for (rid=0; rid < contigs.size(); rid++) {
            if (ref == contigs[rid].first) {
                break;
            }
        }
        if (rid == contigs.size()) {
            return Status::Invalid("range_of_yaml: unknown contig", ref);
        }
    } else if (default_rid >= 0) {
        rid = default_rid;
    } else {
        return Status::Invalid("range_of_yaml: missing 'ref' field");
    }

    const auto n_beg = yaml["beg"];
    V(n_beg && n_beg.IsScalar(), "missing/invalid 'beg' field");
    int beg = n_beg.as<int>()-1;

    const auto n_end = yaml["end"];
    V(n_end && n_end.IsScalar(), "missing/invalid 'end' field");
    int end = n_end.as<int>();

    V(beg >= 0 && end >= beg, "invalid beg/end coordinates");

    ans = range(rid, beg, end);
    return Status::OK();
    #undef V
}

Status yaml_of_discovered_alleles(const discovered_alleles& dal,
                                  const std::vector<std::pair<std::string,size_t> >& contigs,
                                  YAML::Emitter& out) {
    Status s;

    out << YAML::BeginSeq;
    for (const auto& p : dal) {
        out << YAML::BeginMap;

        out << YAML::Key << "range" << YAML::Value;
        S(range_yaml(contigs, p.first.pos, out));
        out << YAML::Key << "dna"
            << YAML::Value << p.first.dna;
        out << YAML::Key << "is_ref"
            << YAML::Value << p.second.is_ref;
        out << YAML::Key << "copy_number"
            << YAML::Value << p.second.copy_number;

        out << YAML::EndMap;
    }
    out << YAML::EndSeq;

    return Status::OK();
}

Status discovered_alleles_of_yaml(const YAML::Node& yaml,
                                  const std::vector<std::pair<std::string,size_t> >& contigs,
                                  discovered_alleles& ans) {
    Status s;
    #define V(pred,msg) if (!(pred)) return Status::Invalid("discovered_alleles_of_yaml: " msg)

    V(yaml.IsSequence(), "not a sequence at top level");
    ans.clear();
    for (YAML::const_iterator p = yaml.begin(); p != yaml.end(); ++p) {
        V(p->IsMap(), "invalid entry");

        range rng(-1,-1,-1);
        const auto n_range = (*p)["range"];
        V(n_range, "missing 'range' field in entry");
        S(range_of_yaml(n_range, contigs, rng));
        #define VR(pred,msg) if (!(pred)) return Status::Invalid("discovered_alleles_of_yaml: " msg, rng.str(contigs))

        const auto n_dna = (*p)["dna"];
        VR(n_dna && n_dna.IsScalar(), "missing/invalid 'dna' field in entry");
        const string& dna = n_dna.Scalar();
        VR(dna.size() > 0, "empty 'dna' in entry");
        VR(regex_match(dna, regex_dna), "invalid allele DNA");
        allele al(rng, dna);

        discovered_allele_info ai;
        const auto n_is_ref = (*p)["is_ref"];
        VR(n_is_ref && n_is_ref.IsScalar(), "missing/invalid 'is_ref' field in entry");
        ai.is_ref = n_is_ref.as<bool>();

        const auto n_copy_number = (*p)["copy_number"];
        VR(n_copy_number && n_copy_number.IsScalar(), "missing/invalid 'copy_number' field in entry");
        ai.copy_number = n_copy_number.as<float>();
        VR(std::isfinite(ai.copy_number) && !std::isnan(ai.copy_number), "invalid 'copy_number' field in entry");

        VR(ans.find(al) == ans.end(), "duplicate alleles");
        ans[al] = ai;
        #undef VR
    }
    V(ans.size() >= 1, "not enough alleles");

    #undef V
    return Status::OK();
}

Status unified_site::yaml(const std::vector<std::pair<std::string,size_t> >& contigs,
                          YAML::Emitter& ans) const {
    Status s;

    ans << YAML::BeginMap;

    ans << YAML::Key << "range" << YAML::Value;
    S(range_yaml(contigs, pos, ans));

    if (containing_target.rid >= 0) {
        ans << YAML::Key << "containing_target" << YAML::Value;
        S(range_yaml(contigs, containing_target, ans));
    }

    ans << YAML::Key << "alleles";
    ans << YAML::Value << YAML::Flow << YAML::BeginSeq;
    for (const auto& allele : alleles) {
        ans << allele;
    }
    ans << YAML::EndSeq;

    ans << YAML::Key << "copy_number";
    ans << YAML::Value << YAML::Flow << YAML::BeginSeq;
    for (auto count : copy_number) {
        ans << count;
    }
    ans << YAML::EndSeq;

    ans << YAML::Key << "unification";
    ans << YAML::Value << YAML::BeginSeq;
    // sort unification entries by to, then by range, then by alt
    vector<pair<allele,int>> u(unification.begin(), unification.end());
    sort(u.begin(), u.end(),
         [] (const pair<allele,int>& p1, const pair<allele,int>& p2) {
            if (p1.second != p2.second) return p1.second < p2.second;
            if (p1.first.pos != p2.first.pos) return p1.first.pos < p2.first.pos;
            return p1.first.dna < p2.first.dna;
         });
    for (const auto& p : u) {
        ans << YAML::BeginMap;
        ans << YAML::Key << "range" << YAML::Value;
        S(range_yaml(contigs, p.first.pos, ans, true));
        ans << YAML::Key << "alt";
        ans << YAML::Value << p.first.dna;
        ans << YAML::Key << "to";
        ans << YAML::Value << p.second;
        ans << YAML::EndMap;
    }
    ans << YAML::EndSeq;

    ans << YAML::EndMap;
    return Status::OK();
}

Status unified_site::of_yaml(const YAML::Node& yaml, const vector<pair<string,size_t> >& contigs,
                             unified_site& ans) {
    Status s;
    #define V(pred,msg) if (!(pred)) return Status::Invalid("unified_site_of_yaml: " msg)

    V(yaml.IsMap(), "not a map at top level");

    const auto n_range = yaml["range"];
    V(n_range, "missing 'range' field");
    S(range_of_yaml(n_range, contigs, ans.pos));
    #define VR(pred,msg) if (!(pred)) return Status::Invalid("unified_site_of_yaml: " msg, ans.pos.str(contigs))

    const auto n_containing_target = yaml["containing_target"];
    if (n_containing_target) {
        S(range_of_yaml(n_containing_target, contigs, ans.containing_target));
    }

    ans.alleles.clear();
    const auto n_alleles = yaml["alleles"];
    VR(n_alleles && n_alleles.IsSequence(), "missing 'alleles' field");
    for (YAML::const_iterator al = n_alleles.begin(); al != n_alleles.end(); ++al) {
        V(al->IsScalar(), "invalid allele");
        ans.alleles.push_back(al->Scalar());
    }
    VR(ans.alleles.size() >= 2, "not enough alleles");

    ans.unification.clear();
    const auto n_unification = yaml["unification"];
    VR(n_unification && n_unification.IsSequence(), "missing 'unification' field");
    for (YAML::const_iterator p = n_unification.begin(); p != n_unification.end(); ++p) {
        VR(p->IsMap(), "invalid unification entry");

        range urange(-1,-1,-1);
        const auto n_urange = (*p)["range"];
        VR(n_urange, "missing 'range' field in unification entry");
        S(range_of_yaml(n_urange, contigs, urange, ans.pos.rid));
        VR(ans.pos.overlaps(urange), "unification entry range does not overlap site range");

        const auto n_ualt = (*p)["alt"];
        VR(n_ualt && n_ualt.IsScalar(), "missing/invalid 'alt' field in unification entry");
        const string& alt = n_ualt.Scalar();
        VR(alt.size() > 0, "empty 'alt' in unification entry");

        const auto n_uto = (*p)["to"];
        VR(n_uto && n_uto.IsScalar(), "missing/invalid 'to' field in unification entry");
        int to = n_uto.as<int>();
        VR(to >= 0 && to < ans.alleles.size(), "invalid 'to' field in unification entry");

        allele al(urange,alt);
        VR(ans.unification.find(al) == ans.unification.end(), "duplicate unification entries");
        ans.unification[al] = to;
    }
    VR(ans.unification.size() >= 2, "not enough unification entries");

    ans.copy_number.clear();
    const auto n_obs = yaml["copy_number"];
    VR(n_obs && n_obs.IsSequence(), "missing 'copy_number' field");
    for (YAML::const_iterator ct = n_obs.begin(); ct != n_obs.end(); ++ct) {
        VR(ct->IsScalar(), "invalid copy_nuber");
        float ctf = ct->as<float>();
        VR(ctf == ctf && ctf >= 0.0, "invalid copy_number");
        ans.copy_number.push_back(ctf);
    }
    VR(ans.copy_number.size() == ans.alleles.size(), "copy_number list has wrong length");

    #undef V
    #undef VR
    return Status::OK();
}

Status unifier_config::of_yaml(const YAML::Node& yaml, unifier_config& ans) {
    Status s;

    ans = unifier_config();

    #define V(pred,msg) if (!(pred)) return Status::Invalid("unifier_config::of_yaml: " msg)
    V(yaml.IsMap(), "not a map at top level");

    const auto n_max_alleles_per_site = yaml["max_alleles_per_site"];
    if (n_max_alleles_per_site) {
        V(n_max_alleles_per_site.IsScalar(), "invalid max_alleles_per_site");
        int max_alleles_per_site = n_max_alleles_per_site.as<int>();
        V(max_alleles_per_site > 0, "invalid max_alleles_per_site");
        ans.max_alleles_per_site = (size_t) max_alleles_per_site;
    }

    const auto n_preference = yaml["preference"];
    if (n_preference) {
        V(n_preference.IsScalar(), "invalid preference");
        const string& preference = n_preference.Scalar();
        if (preference == "common") {
            ans.preference = UnifierPreference::Common;
        } else if (preference == "small") {
            ans.preference = UnifierPreference::Small;
        } else {
            V(false, "invalid preference");
        }
    }

    #undef V
    return Status::OK();
}

Status retained_format_field::of_yaml(const YAML::Node& yaml, unique_ptr<retained_format_field>& ans) {
    Status s;
    #define V(pred,msg) if (!(pred)) return Status::Invalid("retained_format_field::of_yaml: " msg);

    const auto n_orig_names = yaml["orig_names"];
    V(n_orig_names && (n_orig_names.size() > 0), "missing orig_names");
    vector<string> orig_names;
    for (YAML::const_iterator it = n_orig_names.begin(); it != n_orig_names.end(); ++it) {
        V(it->IsScalar(), "invalid orig_names");
        orig_names.push_back(it->Scalar());
    }

    string name;
    const auto n_name = yaml["name"];
    V(n_name && (n_name.Scalar().size() > 0), "missing name");
    name = n_name.Scalar();

    string description;
    const auto n_description = yaml["description"];
    V(n_description && (n_description.Scalar().size() > 0), "missing description");
    description = n_description.Scalar();

    RetainedFieldType type;
    const auto n_type = yaml["type"];
    V(n_type && n_type.IsScalar(), "missing type");
    string s_type = n_type.Scalar();
    if (s_type == "int") {
        type = RetainedFieldType::INT;
    } else if (s_type == "float") {
        type = RetainedFieldType::FLOAT;
    } else {
        V(false, "invalid type");
    }

    RetainedFieldNumber number;
    const auto n_number = yaml["number"];
    V(n_number && n_number.IsScalar(), "missing number");
    string s_number = n_number.Scalar();
    if (s_number == "basic") {
        number = RetainedFieldNumber::BASIC;
    } else if (s_number == "alt") {
        number = RetainedFieldNumber::ALT;
    } else if (s_number == "genotype") {
        number = RetainedFieldNumber::GENOTYPE;
    } else if (s_number == "alleles") {
        number = RetainedFieldNumber::ALLELES;
    } else {
        V(false, "invalid number");
    }

    int count = 0;
    if (number == RetainedFieldNumber::BASIC) {
        const auto n_count = yaml["count"];
        V(n_count && n_count.IsScalar(), "missing count");
        count = n_count.as<int>();
    }

    bool default_to_zero = false;
    const auto n_default_to_zero = yaml["default_to_zero"];
    if (n_default_to_zero) {
        V(n_default_to_zero.IsScalar(), "invalid default_to_zero value");
        string s_default_to_zero = n_default_to_zero.Scalar();

        if (s_default_to_zero == "true") {
            default_to_zero = true;
        } else if (s_default_to_zero == "false") {
            default_to_zero = false;
        } else {
            V(false, "invalid default_to_zero value")
        }
    }

    FieldCombinationMethod combi_method;
    const auto n_combi_method = yaml["combi_method"];
    V(n_combi_method && n_combi_method.IsScalar(), "missing combi_method");
    string s_combi_method = n_combi_method.Scalar();
    if (s_combi_method == "min") {
        combi_method = FieldCombinationMethod::MIN;
    } else if (s_combi_method == "max") {
        combi_method = FieldCombinationMethod::MAX;
    } else {
        V(false, "invalid combi_method");
    }

    ans.reset(new retained_format_field(orig_names, name, type, combi_method, number, count, default_to_zero));
    ans->description = description;
    #undef V
    return Status::OK();
}

Status genotyper_config::of_yaml(const YAML::Node& yaml, genotyper_config& ans) {
    Status s;
    ans = genotyper_config();
    #define V(pred,msg) if (!(pred)) return Status::Invalid("genotyper_config::of_yaml: " msg);

    const auto n_required_dp = yaml["required_dp"];
    if (n_required_dp) {
        V(n_required_dp.IsScalar(), "invalid required_dp");
        ans.required_dp = n_required_dp.as<int>();
    }

    // TODO: allele_dp_format, ref_symbolic_allele, ref_dp_format
    // TODO: output_residuals, output_format

    const auto n_liftover_fields = yaml["liftover_fields"];
    if (n_liftover_fields) {
        V(n_liftover_fields.IsSequence(), "invalid liftover_fields");
        for (YAML::const_iterator it = n_liftover_fields.begin(); it != n_liftover_fields.end(); ++it) {
            unique_ptr<retained_format_field> rff;
            S(retained_format_field::of_yaml(*it, rff));
            ans.liftover_fields.push_back(*rff);
        }
    }

    #undef V
    return Status::OK();
}

} // namespace GLnexus
