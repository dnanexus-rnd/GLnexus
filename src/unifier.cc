#include <algorithm>
#include <assert.h>
#include "unifier.h"
#include <iostream>

using namespace std;

namespace GLnexus {

using discovered_allele = pair<allele,discovered_allele_info>;

// ALT alleles discovered across numerous samples may represent the same edit
// to the reference genome in different ways, specifically if they have
// different amounts of reference padding on either end. (Individual gVCF ALT
// alleles may be reference-padded to unify their representation with other
// ALT alleles observed in the same sample.) So we introduce a notion of
// "minimized" ALT allele, which strips any reference padding so that
// equivalent ALT alleles can be combined, while still remembering the
// original representations.
struct minimized_allele_info {
    set<allele> originals;
    float observation_count = 0.0;

    string str() const {
        ostringstream os;
        os << "Minimized from originals: " << endl;
        for (auto& al : originals) {
            os << "  " << al.str() << endl;
        }
        os << "Observation count: " << observation_count << endl;
        return os.str();
    }
};
using minimized_alleles = map<allele,minimized_allele_info>;
using minimized_allele = pair<allele,minimized_allele_info>;

/// Minimize an ALT allele by trimming bases equal to the reference from
/// either end. This assumes indel alleles are already left-aligned.
Status minimize_allele(const allele& ref, allele& alt) {
    if (ref.pos != alt.pos) return Status::Invalid("minimize_allele: wrong ref allele provided");

    size_t l,r;

    // trim from the right first, to maintain indel left-alignment
    for (r=0; r<alt.dna.size()-1 && r<ref.dna.size()-1
              && alt.dna[alt.dna.size()-r-1] == ref.dna[ref.dna.size()-r-1]; r++);

    // then trim from the left if possible
    for (l=0; l<alt.dna.size()-r-1 && l<ref.dna.size()-r-1 && alt.dna[l] == ref.dna[l]; l++);
    assert(l+r < alt.dna.size());
    assert(l+r < ref.dna.size());

    alt.pos.beg += l;
    alt.pos.end -= r;
    assert(alt.pos.beg <= alt.pos.end);
    alt.dna = alt.dna.substr(l, alt.dna.size()-l-r);
    assert(alt.dna.size());

    return Status::OK();
}

// Separate discovered alleles into the REF alleles and minimized ALT alleles
Status minimize_alleles(const discovered_alleles& src,
                        discovered_alleles& refs, minimized_alleles& alts) {
    Status s;

    // separate the ref and alt alleles
    map<range,discovered_allele> refs_by_range;
    discovered_alleles dalts;
    refs.clear();
    for (const auto& dal : src) {
        if (dal.second.is_ref) {
            refs.insert(dal);
            refs_by_range.insert(make_pair(dal.first.pos, dal));
        } else {
            dalts.insert(dal);
        }
    }

    // go through each alt allele
    alts.clear();
    for (const auto& dal : dalts) {
        allele alt = dal.first;

        // find corresponding reference allele
        const auto rp = refs_by_range.find(alt.pos);
        if (rp == refs_by_range.end()) return Status::Invalid("minimize_alleles: missing REF allele for ", alt.pos.str());

        // minimize the alt allele
        S(minimize_allele(rp->second.first, alt));

        // add it to alts, combining originals and observation_count with any
        // previously observed occurrences of the same minimized alt allele.
        auto ap = alts.find(alt);
        if (ap == alts.end()) {
            minimized_allele_info info;
            info.originals.insert(dal.first);
            info.observation_count = dal.second.observation_count;
            alts[alt] = move(info);
        } else {
            ap->second.originals.insert(dal.first);
            ap->second.observation_count += dal.second.observation_count;
        }
    }

    return Status::OK();
}


// Partition alleles into non-overlapping "active regions" where all the
// alleles in an active region are transitively connected through overlap.
// (However, individual pairs of alleles within an active region might be non-
// overlapping.)
template<class discovered_or_minimized_alleles>
auto partition(const discovered_or_minimized_alleles& alleles) {
    map<range,discovered_or_minimized_alleles> ans;

    range rng(-1,-1,-1);
    discovered_or_minimized_alleles als;

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

// Given a cluster of related/overlapping ALT alleles, prune some so that all
// remaining alleles further cluster into "sites" in which all alleles share
// at least one reference position. We heuristically try to avoid pruning the
// most common alleles. (We have an 0-1 knapsack problem where alleles are
// valued by frequency, and alleles not sharing at least one reference
// position in common conflict.)
auto prune_alleles(const minimized_alleles& alleles, minimized_alleles& pruned) {
    // Sort the alleles by decreasing observation_count.
    vector<minimized_allele> valleles(alleles.begin(), alleles.end());
    sort(valleles.begin(), valleles.end(),
         [] (const minimized_allele& p1, const minimized_allele& p2) {
           return p2.second.observation_count < p1.second.observation_count;
         });

    // Build the sites by considering each alt allele in that order, and
    // accepting or rejecting it as follows. If the allele overlaps exactly
    // one existing site, and overlaps all existing alleles at that site, then
    // accept it. Also accept the allele if it doesn't overlap any already-
    // accepted allele, creating a new site. Reject alleles that overlap some
    // but not all alleles at a site, or that overlap multiple existing sites.
    pruned.clear();
    map<range,minimized_alleles> sites;
    for (const minimized_allele& mal : valleles) {
        // find existing site(s) overlapping this allele
        auto related_site = sites.end();
        bool multiple_sites = false;
        for (auto site = sites.begin(); site != sites.end(); site++) {
            if (mal.first.pos.overlaps(site->first)) {
                if (related_site == sites.end()) {
                    related_site = site;
                } else {
                    multiple_sites = true;
                    break;
                }
            }
        }
        if (multiple_sites) {
            // reject since this allele overlaps multiple existing sites
            pruned.insert(mal);
            continue;
        }
        if (related_site != sites.end()) {
            if (any_of(related_site->second.begin(), related_site->second.end(),
                       [&mal] (const minimized_allele& other) {
                            return !(mal.first.pos.overlaps(other.first.pos));
                        })) {
                // reject since this allele doesn't overlap all alleles in the site
                pruned.insert(mal);
                continue;
            }

            // merge this allele into the related site
            range urng(mal.first.pos.rid,
                       min(mal.first.pos.beg, related_site->first.beg),
                       max(mal.first.pos.end, related_site->first.end));
            minimized_alleles mals(move(related_site->second));
            mals.insert(mal);
            sites.erase(related_site);
            sites[urng] = move(mals);
        } else {
            // no related site: insert a new site
            sites[mal.first.pos] = minimized_alleles({mal});
        }
    }

    return sites;
}

// Given all discovered alleles, delineate one or more sites, where all
// alleles at a site share at least one position; some alleles may be pruned
// to achieve this. The result maps the range of each site to a pair of the
// discovered REF alleles and the minimized ALT alleles.
Status delineate_sites(const discovered_alleles& alleles,
                       map<range,pair<discovered_alleles,minimized_alleles>>& ans) {
    Status s;

    // Start by finding non-interacting "active regions" -- individual pairs
    // of alleles within each active region might be non-overlapping.
    auto active_regions = partition<discovered_alleles>(alleles);

    // Decompose each active region into sites in which all alleles share at
    // least one position.
    ans.clear();
    for (const auto& active_region : active_regions) {
        // minimize the alt alleles
        discovered_alleles refs;
        minimized_alleles alts, pruned;
        S(minimize_alleles(active_region.second, refs, alts));

        // prune alt alleles as necessary to yield sites
        auto sites = prune_alleles(alts, pruned);
        // TODO: should we do anything with the pruned alleles?

        for (const auto& site : sites) {
            // find the ref alleles overlapping this site
            discovered_alleles site_refs;
            for (const auto& ref : refs) {
                if (ref.first.pos.overlaps(site.first)) {
                    site_refs.insert(ref);
                }
            }
            assert(ans.find(site.first) == ans.end());
            ans[site.first] = make_pair(site_refs,site.second);
        }
    }
    return Status::OK();
}

// Given the alleles at a site, collapse the discovered REF alleles to get one
// reference allele covering the site's range.
Status unify_ref(const range& pos, const discovered_alleles& refs, allele& ref) {
    ref.pos = pos;
    ref.dna.assign(pos.size(), char(0));

    for (const auto& ref1 : refs) {
        UNPAIR(ref1, al, info);
        assert(pos.overlaps(al.pos));
        assert(al.pos.size() == al.dna.size());
        for (int i = 0, j = (al.pos.beg-pos.beg); i < al.dna.size(); i++, j++) {
            if (j >= 0 && j < ref.dna.size()) {
                if (ref.dna[j] && ref.dna[j] != al.dna[i]) {
                    return Status::Invalid("detected inconsistent REF alleles", al.pos.str());
                }
                ref.dna[j] = al.dna[i];
            }
        }
    }
    if (any_of(ref.dna.begin(), ref.dna.end(), [] (const char c) { return c == 0; })) {
        return Status::Invalid("Incomplete REF allele coverage in unification (this should not happen!)", pos.str());
    }

    return Status::OK();
}

// Given the unified REF allele and an ALT allele in the same site, pad the
// ALT allele with reference bases so that it covers exactly the site's range.
Status pad_alt_allele(const allele& ref, allele& alt) {
    if (!ref.pos.contains(alt.pos)) return Status::Invalid("BUG: unified REF allele doesn't contain all ALT alleles");
    size_t l = alt.pos.beg - ref.pos.beg;
    size_t r = ref.pos.end - alt.pos.end;

    alt.pos = ref.pos;
    alt.dna = ref.dna.substr(0, l) + alt.dna + ref.dna.substr(ref.dna.size()-r, r);
    return Status::OK();
}

/// Unify the alleles at one site. The given ALT alleles must all share at
/// least one position in common.
Status unify_alleles(const unifier_config& cfg, const range& pos,
                     const discovered_alleles& refs, const minimized_alleles& alts,
                     unified_site& ans) {
    Status s;

    // collapse the refs to get the reference allele for this site
    allele ref(pos, "A");
    S(unify_ref(pos, refs, ref));

    // sort the alt alleles by decreasing observation count, then by minimized
    // DNA size, then by minimized position, then lexicographically by DNA
    vector<minimized_allele> valts(alts.begin(), alts.end());
    sort(valts.begin(), valts.end(),
         [] (const minimized_allele& p1, const minimized_allele& p2) {
            assert (p1.first.pos.overlaps(p2.first.pos));
            if (p2.second.observation_count != p1.second.observation_count)
                return p2.second.observation_count < p1.second.observation_count;
            if (p1.first.dna.size() != p2.first.dna.size())
                return p1.first.dna.size() < p2.first.dna.size();
            if (p1.first.pos != p2.first.pos)
                return p1.first.pos < p2.first.pos;
            return p1.first.dna < p2.first.dna;
         });

    // enforce max_alleles_per_site. NB, valts does not include the ref
    // allele, so we're truncating valts to cfg.max_alleles_per_site-1
    if (cfg.max_alleles_per_site > 1 && valts.size() >= cfg.max_alleles_per_site) {
        // Out of fairness, we also prune alt alleles with the same
        // observation count as the first truncated one, except we always need
        // to keep at least one alt allele of course.
        unsigned int trunc0 = 1, trunc1 = 1;
        for (; trunc1 < cfg.max_alleles_per_site; trunc1++) {
            assert(valts[trunc0].second.observation_count >= valts[trunc1].second.observation_count);
            while(valts[trunc0].second.observation_count > valts[trunc1].second.observation_count)
                trunc0++;
        }
        // TODO: should we do something with the additional pruned alleles?
        valts.assign(valts.begin(), valts.begin()+trunc0);
        assert(valts.size() < cfg.max_alleles_per_site);
    }

    // fill out the unification
    unified_site us(pos);
    us.alleles.push_back(ref.dna);
    // We don't attempt to specify the ref observation count for now. It's
    // problematic because the observation counts in discovered_alleles omit
    // nearly all homozygous ref genotypes.
    us.observation_count.push_back(0.0);
    for (const auto& ref : refs) {
        us.unification[ref.first] = 0;
    }
    for (int i = 1; i <= valts.size(); i++) {
        const minimized_allele& p = valts[i-1];
        UNPAIR(p, alt, alt_info);

        allele unified_alt(alt);
        S(pad_alt_allele(ref, unified_alt));
        assert(unified_alt.pos == ref.pos);

        us.alleles.push_back(unified_alt.dna);
        us.observation_count.push_back(alt_info.observation_count);
        for (const auto& original : alt_info.originals) {
            us.unification[original] = i;
        }
    }
    ans = std::move(us);
    return Status::OK();
    // TODO: will need some representation for observation count of pruned
    // alleles. but maybe counting samples (+ploidy config) in
    // discovered_alleles will suffice.
}

Status unified_sites(const unifier_config& cfg,
                     const discovered_alleles& alleles,
                     vector<unified_site>& ans) {
    Status s;

    map<range,pair<discovered_alleles,minimized_alleles>> sites;
    S(delineate_sites(alleles, sites));

    ans.clear();
    for (const auto& site : sites) {
        UNPAIR(site, pos, site_alleles);
        UNPAIR(site_alleles, ref_alleles, alt_alleles);
        unified_site us(pos);
        S(unify_alleles(cfg, pos, ref_alleles, alt_alleles, us));
        ans.push_back(us);
    }
    return Status::OK();
}

}
