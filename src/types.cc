#include "types.h"

using namespace std;

namespace GLnexus {

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

Status unified_site_of_yaml(const YAML::Node& yaml, const vector<pair<string,size_t> >& contigs,
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
        VR(ans.pos.contains(urange), "unification entry range isn't contained within site range");

        const auto n_ualt = (*p)["alt"];
        VR(n_ualt && n_ualt.IsScalar(), "missing/invalid 'alt' field in unification entry");
        const string& alt = n_ualt.Scalar();
        VR(alt.size() > 0, "empty 'alt' in unification entry");

        const auto n_uto = (*p)["to"];
        VR(n_uto && n_uto.IsScalar(), "missing/invalid 'to' field in unification entry");
        int to = n_uto.as<int>();
        VR(to >= 0 && to < ans.alleles.size(), "invalid 'to' field in unification entry");

        auto k = make_pair(urange,alt);
        VR(ans.unification.find(k) == ans.unification.end(), "duplicate unification entries");
        ans.unification[k] = to;
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

    return Status::OK();
}

}
