#include <algorithm>
#include <assert.h>
#include "alleles.h"
#include <iostream>

using namespace std;

namespace GLnexus {

using discovered_allele = pair<allele,discovered_allele_info>;

// Partition discovered alleles into non-overlapping "sites". However,
// individual pairs of alleles within "sites" might be non-overlapping
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
// alleles at a site share at least one position. Alleles may be pruned to
// satisfy this. We can view this as an 0-1 knapsack problem with values and
// conflicts. We heuristically try to avoid pruning the most common alleles.
auto prune_alleles(const discovered_alleles& alleles, discovered_alleles& pruned) {
    map<range,discovered_alleles> ans;
    pruned.clear();

    // First coarsely partition the alt alleles into clusters that don't
    // interact with each other at all
    auto clusters = partition_alleles(alleles);

    // Now decompose each cluster into sites, where all alleles at a site
    // share at least one reference position.
    for (const pair<range,discovered_alleles>& cluster : clusters) {
        // separate the ref and alt alleles
        vector<discovered_allele> ref_alleles, alt_alleles;
        for (const auto& dsal : cluster.second) {
            if (dsal.second.is_ref) {
                ref_alleles.push_back(dsal);
            } else {
                alt_alleles.push_back(dsal);
            }
        }

        // sort alt alleles in the cluster by decreasing observation count
        sort(alt_alleles.begin(), alt_alleles.end(),
             [] (const discovered_allele& p1, const discovered_allele& p2) {
                assert (!p1.second.is_ref && !p2.second.is_ref);
                return p2.second.observation_count < p1.second.observation_count;
             });

        map<range,discovered_alleles> cluster_sites;
        // Build the sites by considering each alt allele and accepting or
        // rejecting it as follows. If the allele overlaps exactly one
        // existing site, accept it if it overlaps all existing alleles at
        // that site. Also accept the allele if it doesn't overlap any
        // already-accepted allele, creating a new site. Reject alleles that
        // overlap some but not all alleles at a site, or that overlap
        // multiple existing sites.
        for (const discovered_allele& dal : alt_alleles) {
            // find existing site(s) overlapping this allele
            auto related_site = cluster_sites.end();
            bool multiple_sites = false;
            for (auto site = cluster_sites.begin(); site != cluster_sites.end(); site++) {
                if (dal.first.pos.overlaps(site->first)) {
                    if (related_site == cluster_sites.end()) {
                        related_site = site;
                    } else {
                        multiple_sites = true;
                    }
                }
            }
            if (multiple_sites) {
                // reject since this allele overlaps multiple existing sites
                pruned.insert(dal);
                continue;
            }
            if (related_site != cluster_sites.end()) {
                if (any_of(related_site->second.begin(), related_site->second.end(),
                           [&dal] (const discovered_allele& other) {
                                return !(dal.first.pos.overlaps(other.first.pos));
                            })) {
                    // reject since this allele doesn't overlap all alleles in the site
                    pruned.insert(dal);
                    continue;
                }

                // merge this allele into the related site
                range urng(dal.first.pos.rid,
                           min(dal.first.pos.beg, related_site->first.beg),
                           max(dal.first.pos.end, related_site->first.end));
                discovered_alleles dals(related_site->second);
                dals.insert(dal);
                cluster_sites.erase(related_site);
                cluster_sites[urng] = dals;
            } else {
                // no related site: insert a new site
                cluster_sites[dal.first.pos] = discovered_alleles({dal});
            }
        }

        // add back in the ref alleles matching the remaining alt alleles
        for (const auto& refal : ref_alleles) {
            auto site =
                find_if(cluster_sites.begin(), cluster_sites.end(),
                        [&] (const pair<range,discovered_alleles>& p) {
                            return find_if(p.second.begin(), p.second.end(),
                                           [&] (const discovered_allele& dsal) {
                                            return (dsal.first.pos == refal.first.pos);
                                           }) != p.second.end();
                        });
            if (site != cluster_sites.end()) {
                assert(site->first.overlaps(refal.first.pos));
                site->second.insert(refal);
            } else {
                pruned.insert(refal);
            }
        }

        // add the sites from this cluster to the final answer
        for (const auto& site : cluster_sites) {
            ans.insert(site);
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
            if (p1.second.is_ref == true && p2.second.is_ref == false) {
                return true;
            } else if (p2.second.is_ref == true && p1.second.is_ref == false) {
                return false;
            }
            return p2.second.observation_count < p1.second.observation_count;
         });

    // find the most common alt allele
    auto most_common_alt_it =
        find_if(allelesv.begin(), allelesv.end(),
                [&] (const discovered_allele& p) {
                    return p.second.is_ref == false;
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
                    [] (const discovered_allele& p) { return p.second.is_ref == true; }) == 1);
    assert(count_if(allelesv.begin(), allelesv.end(),
                    [] (const discovered_allele& p) { return p.second.is_ref == false; }) > 0);

    // fill out the unification
    unified_site us(most_common_alt.pos);
    for (int i = 0; i < allelesv.size(); i++) {
        const discovered_allele& p = allelesv[i];
        const allele& al = p.first;

        us.alleles.push_back(al.dna);
        us.unification[make_pair(al.pos.beg,al.dna)] = i;
        us.observation_count.push_back(p.second.observation_count);
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
