#include <algorithm>
#include <assert.h>
#include "alleles.h"
#include <iostream>

using namespace std;

namespace GLnexus {

// partition discovered alleles into non-overlapping "sites"
auto partition_alleles(const discovered_alleles& alleles) {
    map<range,discovered_alleles> ans;

    range rng(-1,-1,-1);
    discovered_alleles als;

    for (const auto& it : alleles) {
        if (!rng.overlaps(it.first.pos)) {
            if (rng.rid != -1) {
                assert(!als.empty());
                ans[rng] = als;
            }
            rng = it.first.pos;
            als.clear();
        }

        assert(rng.overlaps(it.first.pos));
        assert(rng <= it.first.pos);
        rng.end = max(rng.end, it.first.pos.end);
        assert(it.first.pos.within(rng));
        als.insert(it);
    }

    if (rng.rid != -1) {
        assert(!als.empty());
        ans[rng] = als;
    }

    return ans;
}

// prune rare alt alleles until the remaining set can be partitioned into non-
// overlapping "sites", where all alleles at a site overlap by at least 1
// position
auto prune_alleles(const discovered_alleles& alleles, discovered_alleles& pruned) {
    map<range,discovered_alleles> ans;
    pruned.clear();

    // separate the ref and alt alleles
    discovered_alleles ref_alleles, working_set;
    for (const auto& dsal : alleles) {
        UNPAIR(dsal,allele,ai)
        bool is_ref = ai.first;
        if (is_ref) {
            ref_alleles.insert(dsal);
        } else {
            working_set.insert(dsal);
        }
    }

    // TODO: smarterer algorithm
    do {
        // partition the alt alleles into non-overlapping sites
        ans = partition_alleles(working_set);

        // at each site, do all alleles overlap by at least 1 position?
        bool solution = true;
        for (const auto& site : ans) {
            UNPAIR(site,site_range,site_alleles)
            for (const auto& dsal : site_alleles) {
                auto common_range = site_range.intersect(dsal.first.pos);
                if (common_range) {
                    site_range = *common_range;
                } else {
                    solution = false;
                    break;
                }
            }
            if (!solution) break;
        }

        // if so, we have a solution
        if (solution) break;

        // otherwise, prune the rarest allele in the working set, and repeat
        assert(working_set.size() > 2);
        auto rare_allele = working_set.cend();
        float rare_count = numeric_limits<float>::max();
        for (auto it = working_set.cbegin(); it != working_set.cend(); it++) {
            if(it->second.second < rare_count) {
                rare_allele = it;
                rare_count = it->second.second;
            }
        }
        assert(rare_allele != working_set.cend());
        pruned.insert(*rare_allele);
        working_set.erase(rare_allele);
    } while(true);

    // add back in the ref alleles matching the remaining alt alleles
    for (const auto& refal : ref_alleles) {
        if (ans.find(refal.first.pos) != ans.end()) {
            ans[refal.first.pos].insert(refal);
        } else {
            pruned.insert(refal);
        }
    }

    return ans;
}

// PLACEHOLDER no-op unification. just groups exact range matches
Status unify_alleles_placeholder(const discovered_alleles& alleles, vector<unified_site>& ans) {
    set<range> ph_ranges;
    multimap<range,tuple<allele,bool,float> > ph_alleles;

    // group alleles by exact range
    for (const auto& it : alleles) {
        UNPAIR(it,allele,ai)
        UNPAIR(ai,is_ref,obs_count)
        ph_ranges.insert(allele.pos);
        ph_alleles.insert(make_pair(allele.pos,make_tuple(allele,is_ref,obs_count)));
    }

    // for each range
    ans.clear();
    for (const auto& ph_range : ph_ranges) {
        unified_site us(ph_range);

        // fill in the list of allele DNAs, with the ref allele first
        string ph_ref_dna;
        float ph_ref_obs;
        vector<pair<string,float> > us_alleles;
        auto its = ph_alleles.equal_range(ph_range);
        for (auto it = its.first; it != its.second; it++) {
            UNPAIR(*it,range2,p)
            assert(ph_range == range2);
            auto allele = get<0>(p);
            auto is_ref = get<1>(p);
            auto obs_count = get<2>(p);
            if (is_ref) {
                assert(ph_ref_dna.empty());
                ph_ref_dna = allele.dna;
                ph_ref_obs = obs_count;
            } else {
                us_alleles.push_back(make_pair(allele.dna,obs_count));
            }
        }
        assert(!ph_ref_dna.empty());
        us.alleles.push_back(ph_ref_dna);
        us.observation_count.push_back(ph_ref_obs);

        // sort alt alleles by decreasing observation count
        sort(us_alleles.begin(), us_alleles.end(),
             [] (const pair<string,float>& p1, const pair<string,float>& p2) {
                return get<1>(p2) < get<1>(p1);
             });
        for (const auto& it : us_alleles) {
            us.alleles.push_back(get<0>(it));
            us.observation_count.push_back(get<1>(it));
        }
        

        // fill in the unification
        map<string,int> ph_alleles_indices;
        for (unsigned i = 0; i < us.alleles.size(); i++) {
            ph_alleles_indices[us.alleles[i]] = i;
        }
        for (auto it = its.first; it != its.second; it++) {
            UNPAIR(*it,range2,p)
            assert(ph_range == range2);
            auto allele = get<0>(p);
            us.unification[make_pair(ph_range.beg,allele.dna)] = ph_alleles_indices[allele.dna];
        }

        ans.push_back(us);
    }

    return Status::OK();
}

Status unify_alleles(const discovered_alleles& alleles, vector<unified_site>& ans) {
    return unify_alleles_placeholder(alleles, ans);

    // TODO prune alleles
    // 1. overhanging target region
    // 2. rare
    // 3. top N

    auto sites = partition_alleles(alleles);

    ans.clear();
    for (const auto& site : sites) {
        unified_site us(site.first);

        // compute the reference allele covering the site. all reference sub-alleles unify into it.
        // compute alt alleles by splicing in the non-ref alleles at appropriate offsets; remember their unification.

        for (const auto& allele : site.second) {
            assert(allele.first.pos.within(us.pos));
            us.alleles.push_back(allele.first.dna);
            // TODO: actually unify alleles!
            // TODO: ensure ref allele comes first
        }
        ans.push_back(us);
    }

    return Status::OK();
}

}
