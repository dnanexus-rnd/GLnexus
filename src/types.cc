#include "types.h"
#include <algorithm>

using namespace std;

namespace GLnexus {

bool is_dna(const string& str) {
    return all_of(str.begin(), str.end(),
                  [](char ch) { return ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch =='N'; });
}

// Add src alleles to dest alleles. Identical alleles alleles are merged,
// using the sum of their observation_counts
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
            p->second.observation_count += ai.observation_count;
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

    V(beg >= 1 && end >= beg, "invalid beg/end coordinates");

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
        out << YAML::Key << "observation_count"
            << YAML::Value << p.second.observation_count;

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
        VR(is_dna(dna), "invalid allele DNA");
        allele al(rng, dna);

        discovered_allele_info ai;
        const auto n_is_ref = (*p)["is_ref"];
        VR(n_is_ref && n_is_ref.IsScalar(), "missing/invalid 'is_ref' field in entry");
        ai.is_ref = n_is_ref.as<bool>();

        const auto n_observation_count = (*p)["observation_count"];
        VR(n_observation_count && n_observation_count.IsScalar(), "missing/invalid 'observation_count' field in entry");
        ai.observation_count = n_observation_count.as<float>();
        VR(std::isfinite(ai.observation_count) && !std::isnan(ai.observation_count), "invalid 'observation_count' field in entry");

        VR(ans.find(al) == ans.end(), "duplicate alleles");
        ans[al] = ai;
        #undef VR
    }
    V(ans.size() >= 1, "not enough alleles");

    #undef V
    return Status::OK();
}

Status unified_site::yaml(const std::vector<std::pair<std::string,size_t> >& contigs,
                          YAML::Emitter& ans) {
    Status s;

    ans << YAML::BeginMap;

    ans << YAML::Key << "range" << YAML::Value;
    S(range_yaml(contigs, pos, ans));

    ans << YAML::Key << "alleles";
    ans << YAML::Value << YAML::Flow << YAML::BeginSeq;
    for (const auto& allele : alleles) {
        ans << allele;
    }
    ans << YAML::EndSeq;

    ans << YAML::Key << "observation_count";
    ans << YAML::Value << YAML::Flow << YAML::BeginSeq;
    for (auto count : observation_count) {
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

    ans.observation_count.clear();
    const auto n_obs = yaml["observation_count"];
    VR(n_obs && n_obs.IsSequence(), "missing 'observation_count' field");
    for (YAML::const_iterator ct = n_obs.begin(); ct != n_obs.end(); ++ct) {
        VR(ct->IsScalar(), "invalid observation count");
        float ctf = ct->as<float>();
        VR(ctf == ctf && ctf >= 0.0, "invalid observation_count");
        ans.observation_count.push_back(ctf);
    }
    VR(ans.observation_count.size() == ans.alleles.size(), "observation_count list has wrong length");

    #undef V
    #undef VR
    return Status::OK();
}

}
