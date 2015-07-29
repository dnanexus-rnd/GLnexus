#include <algorithm>
#include <assert.h>
#include "alleles.h"
#include <iostream>

using namespace std;

namespace GLnexus {

using discovered_allele = pair<allele,pair<bool,float>>;

// Partition discovered alleles into non-overlapping "sites". However,
// individual pairs of alleles within "sites" may not overlap.
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

// Partition discovered alleles into non-overlapping "sites", where all
// alleles at a site share at least one position. Prunes rare alleles until
// that condition can be satisfied.
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
        auto site =
            find_if(ans.begin(), ans.end(),
                    [&] (const pair<range,discovered_alleles>& p) {
                        return find_if(p.second.begin(), p.second.end(),
                                       [&] (const discovered_allele& dsal) {
                                        return (dsal.first.pos == refal.first.pos);
                                       }) != p.second.end();
                    });
        if (site != ans.end()) {
            assert(site->first.overlaps(refal.first.pos));
            site->second.insert(refal);
        } else {
            pruned.insert(refal);
        }
    }

    return ans;
}

Status delineate_sites(const discovered_alleles& alleles, map<range,discovered_alleles>& ans) {
    // Start with an initial partition_alleles to delineate totally non-
    // interacting active regions. Individual pairs of alleles within each
    // active region may be non-overlapping.
    auto active_regions = partition_alleles(alleles);

    // Use prune_alleles to decompose each active region into sites in which
    // all alleles share at least one position -- rare or lengthy alleles may
    // need to be pruned in order to achieve this
    ans.clear();
    for (const auto& active_region : active_regions) {
        discovered_alleles pruned_alleles;
        auto site = prune_alleles(active_region.second, pruned_alleles);
        ans.insert(site.begin(), site.end());
    }
    return Status::OK();
}

// Placeholder unification: avoids the need for reference padding by pruning
// any alleles covering a different reference range than the most common
// alt allele
//
// alleles is a set of alleles which all overlap by at least 1 position
Status unify_alleles_placeholder(const discovered_alleles& alleles, unified_site& ans) {

    // sort the alleles ref first, then by decreasing observation count
    vector<discovered_allele> allelesv(alleles.begin(), alleles.end());
    sort(allelesv.begin(), allelesv.end(),
         [] (const discovered_allele& p1, const discovered_allele& p2) {
            if (p1.second.first == true && p2.second.first == false) {
                return true;
            } else if (p2.second.first == true && p1.second.first == false) {
                return false;
            }
            return p2.second.second < p1.second.second;
         });

    // find the most common alt allele
    auto most_common_alt_it =
        find_if(allelesv.begin(), allelesv.end(),
                [&] (const discovered_allele& p) {
                    return p.second.first == false;
                });
    assert(most_common_alt_it != allelesv.end());
    auto most_common_alt = most_common_alt_it->first;

    // remove alleles not covering exactly the same reference range
    allelesv.erase(remove_if(allelesv.begin(), allelesv.end(),
                             [&] (const discovered_allele& p) {
                                return p.first.pos != most_common_alt.pos;
                             }), allelesv.end());

    // sanity check: exactly one reference allele should remain
    assert(count_if(allelesv.begin(), allelesv.end(),
                    [] (const discovered_allele& p) { return p.second.first == true; }) == 1);
    assert(count_if(allelesv.begin(), allelesv.end(),
                    [] (const discovered_allele& p) { return p.second.first == false; }) > 0);

    // fill out the unification
    unified_site us(most_common_alt.pos);
    for (int i = 0; i < allelesv.size(); i++) {
        const discovered_allele& p = allelesv[i];
        const allele& al = p.first;

        us.alleles.push_back(al.dna);
        us.unification[make_pair(al.pos.beg,al.dna)] = i;
        us.observation_count.push_back(p.second.second);
    }
    ans = std::move(us);
    return Status::OK();
}

/// Unify the alleles at one site
Status unify_alleles(const range& pos, const discovered_alleles& alleles, unified_site& ans) {
    return unify_alleles_placeholder(alleles, ans);
}

Status unified_sites(const discovered_alleles& alleles, vector<unified_site>& ans) {
    Status s;

    map<range,discovered_alleles> sites;
    S(delineate_sites(alleles, sites));

    ans.clear();
    for (const auto& site : sites) {
        UNPAIR(site, pos, site_alleles);
        unified_site us(pos);
        S(unify_alleles(pos, site_alleles, us));
        ans.push_back(us);
    }
    return Status::OK();
}

}
